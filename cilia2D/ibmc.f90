program ibmc
    use mod_pressure,       only: generate_laplacian_sparse, generate_laplacian_sparse_constant, calculate_pressure_sparse
    use mod_amgx,           only: calculate_pressure_amgx
    use mod_mesh
    use mod_time
    use mod_boundary
    use mod_io
    use mod_ibm
    use fem2d
    use mod_cilia
    use mod_particle
    implicit none

    ! AmgX
    logical         :: init_status 

    ! Parameters
    real(8), parameter :: PI = 3.141592653589793
    ! Computational Domain
    real(8)    :: Lx
    real(8)    :: Ly
    ! Mesh Paramaters
    integer  :: Nx
    integer  :: Ny
    ! Simulation time Paramaters
    real(8)    :: tsim
    real(8)    :: dt  
    real(8)    :: t
    integer  :: it_save
    ! Physical Constants
    real(8)    :: nu 
    real(8)    :: rho   
    ! Boundary values
    real(8)    :: utop, vtop, ubottom, vbottom, &
                       uleft, vleft, uright, vright, ulid
    ! Matrices to store fields
    real(8), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), &
                                 P(:,:), A(:,:,:), Fx(:,:), Fy(:,:), umid(:,:), vmid(:,:), w(:,:), &
                                 au(:,:), av(:,:), du(:,:), dv(:,:), Pold(:,:), dpx(:,:), dpy(:,:)
    ! Temporary/Miscellaneous variable
    integer  :: it, NN, ic, err
    integer :: niter, iter
    real(8) :: tp  ! Time period
    real(8) :: tbegin
    ! Mesh
    type(mesh)      :: M
    ! AB: crank nicolson factor
    real(real64) :: cnfac, dxi, dyi, aij

    ! Cilia information
    real(8) :: x0(2)
    real(8) :: ftip  
    real(8) :: mtip  
    integer :: ncilia
    integer :: icilia
    real(8) :: L
    integer :: j, i
    integer :: nvc ! Number of vertices in the cilia
    ! type(cilium), allocatable :: cilia(:)
    type(cilium), allocatable :: cilia(:)
    real(8) :: Kc, Bc, kap, bap
    real(8), allocatable :: xcil(:,:)
    real(8) :: lcbed, delc
    
    ! ! Particle information
    ! integer :: nparticle, num_rows
    ! integer :: nvp ! Number of vertices in the particle
    ! type(particle), allocatable :: particles(:)
    ! real(real64), allocatable :: thetav(:)
    integer :: nvp ! Number of vertices in the particle
    real(real64) :: aa, bb ! Major and minor axis
    real(8) :: Kp, Bp

    ! FEM data
    type(festruct), allocatable :: particles(:)
    integer(int32) :: ntri, npp, nparticle
    real(real64), allocatable :: pb(:,:,:)
    integer, allocatable :: mp(:,:)
    real(real64), allocatable :: paelem(:)
    real(real64), allocatable :: pp(:,:)
    real(real64), allocatable :: FN(:,:)
    real(real64), allocatable :: UN(:,:)
    logical, allocatable :: pboundary(:,:)
    integer(int32) :: femdata(2)
   
    ! Namelists for input
    namelist /time/ dt, it_save, tsim
    namelist /grid/ Nx, Ny, Lx, Ly
    namelist /flow/ nu, rho, utop, TP
    namelist /ciliaprops/ ftip, mtip, ncilia, Kc, nvc, Bc, kap, bap
    namelist /particleprops/ Kp, nvp, Bp, aa, bb

    !---------------------- Begin Calculations ------------------------------------!

    ! Read input data from file
    open(1004,file="input_params.dat",form='formatted')
    READ(unit=1004,nml=grid,iostat=err)
    READ(unit=1004,nml=ciliaprops,iostat=err)
    READ(unit=1004,nml=particleprops,iostat=err)
    close(1004)

    print *, "ftip = ", ftip

    ! Construct and write Mesh data
    M = mesh('M',Lx,Ly,Nx,Ny)
    call write_mesh(M,'u')
    call write_mesh(M,'v')
    call write_mesh(M,'p')

    ! Allocate matrices ofr u,v,us,vs,rhs
    allocate(u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(umid(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(vmid(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(Pold(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(w(M%xo%lb:M%xo%ub,M%yo%lb:M%yo%ub)) ! Vorticity
    NN = 5
    allocate(A(1:Nx,1:Ny,NN))
    allocate(au(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(av(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(du(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(dv(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(dpx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))
    allocate(dpy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))

    ! Allocate fem variables
    !------- femdata----------!
    OPEN(33, FILE="femdata.bin",&
     FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
    READ(33) femdata
    close(33)
    !---------------------!
    npp = femdata(1)
    ntri = femdata(2)
    ! npp = 47!216
    ! ntri = 72!380
    nparticle = 1
    allocate(mp(ntri,3), pb(ntri,2,3), pp(npp,2), paelem(ntri), & 
                FN(npp,2), UN(npp,2), pboundary(npp,2), particles(nparticle))

    print *, npp, ntri

    ! Initialize
    u   = 0.0d0
    v   = 0.0d0
    us  = 0.0d0
    vs  = 0.0d0
    R   = 0.0d0
    A   = 0.0d0
    Fx  = 0.0d0
    Fy  = 0.0d0
    P   = 0.0d0
    w   = 0.0d0
    au  = 0.0d0
    av  = 0.0d0
    du  = 0.0d0
    dv  = 0.0d0
    Pold = 0.0d0


    ! Define boundary conditions for velocity
    utop    = 0.0d0
    vtop    = 0.0d0
    ubottom = 0.0d0
    vbottom = 0.0d0
    uleft   = 0.0d0
    vleft   = 0.0d0
    uright  = 0.0d0
    vright  = 0.0d0

    cnfac=0.0d0

    ! Create cilia
    ! x0 = [Lx/4.0, Ly/8.0]
    x0 = [Lx/4.0, 3*M%dx]
    l = Ly/4.0d0 
    allocate(cilia(ncilia), xcil(ncilia,2))
    lcbed = 0.80d0 * Lx ! Length of cilia bed
    delc = lcbed / (ncilia - 1)
    xcil (:,2) = x0(2)
    xcil(:,1) = (Lx - lcbed)/2.0d0 + delc * [(i-1, i=1, ncilia)]
    ! xcil(:,1) = Lx/2.0
    do i = 1,ncilia
        cilia(i) = cilium(xcil(i,:),l,nvc,Kc,Bc,kap,bap)
    end do
        
    ! Create fem particle
    call read_fem_data(mp,paelem,pb,pp)
    pboundary = .FALSE.
    ! PP(:,1) = PP(:,1) + Lx/2.0d0 
    PP(:,1) = PP(:,1) + xcil(floor(ncilia/2.0),1) ! Place the particle on the top of the middle cilia 
    ! PP(:,2) = PP(:,2) + Ly/2.0d0  
    PP(:,2) = PP(:,2) + (Ly - (3*M%dx + l))/2.0 + (l + 3*M%dx)
    ! kp = kval, co = bp
    particles(1) = festruct(MP,PP,PB,pboundary,paelem,FN,UN,bp,kp,1.0d0) ! kp = kval, co = bp , dl = 1.0d0

    ! Generate Laplacian matrix
    call generate_laplacian_sparse_constant(A,M%dx,M%dy)

    ! Read input data from file
    open(1004,file="input_params.dat",form='formatted')
    READ(unit=1004,nml=time,iostat=err)
    READ(unit=1004,nml=flow,iostat=err)
    close(1004)

    write(*,'(2(A,I8),2(A,1p1e15.6))') "Nx = ",Nx, " Ny = ",Ny, &
    " Lx = ",Lx, " Ly = ",Ly

    !call apply_parabolic_initialization(M,u,utop)
    ! call apply_couette_initialization(M,u,utop*2*3.14/tp,0.0d0*utop)
    print *, ulid
    
    niter=1 ! CN iterations
    
    ! Start time loop
    t = 0.0d0
    it = 0

    init_status = .False. ! AmgX initialization status
    call calculate_rhs(M,u,v,R,rho,dt)
    call calculate_pressure_amgx(A,Pold,R,init_status)

    ! Lid velocity
    ulid = utop
    
    ! to be used for implicit factor
    dxi = 1.0d0/M%dx
    dyi = 1.0d0/M%dy
    
    aij=cnfac*2*(dxi*dxi+dyi*dyi)*nu

    tbegin = 000.0d0!500.0d0
    ! Open the file for writing
    OPEN(UNIT=10, FILE="MP.txt", STATUS='replace', ACTION='write')
    ! Loop through the matrix and write its elements to the file
    DO i = 1, ntri
        WRITE(10, '(3I5)') (particles(1)%M(i, j), j = 1, 3)
    ENDDO
    ! Close the file
    CLOSE(10)

    print *, ftip
    do while (t.lt.tsim)
        ! ulid is the amplitude
        utop = 2*3.1415/tp*ulid*cos(2*3.1415*t/tp)
        ! ubottom = 2*3.1415/tp*ulid*cos(2*3.1415*t/tp)
        ! utop = ulid
        ! ubottom = -ulid
        ! print *, "utop = " , utop
        ! utop = 0.0d0
        t = t + dt
        it = it + 1
      
        write(*,'(A,F18.10)') 'time = ', t

        ! Apply velocity boundary conditions
        call apply_boundary_channel(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

        call particles(:)%calculate_forces()
        call cilia(:)%forces(ftip)
        call cilia(:)%moments()
        Fx = 0.0d0 ! Initialize the forces at every time-step
        Fy = 0.0d0
        call spread_force(M,particles,Fx,Fy)
        call spread_force(M,cilia,Fx,Fy)
        call spread_vorticity_force(M,cilia,Fx,Fy)

        call advection(M,u,v,au,av)
        
        call gradp(M,Pold,dpx,dpy)
        do iter=1,niter
          call diffusion(M,u,v,us,vs,du,dv,cnfac)
        !   us = (u + (nu*du - au + Fx - dpx/rho) * dt)/(1.0d0+aij*dt)
        !   vs = (v + (nu*dv - av + Fy - dpy/rho) * dt)/(1.0d0+aij*dt)

          us = (u + (nu*du - au + Fx) * dt)/(1.0d0+aij*dt)
          vs = (v + (nu*dv - av + Fy) * dt)/(1.0d0+aij*dt)

        ! Apply velocity boundary conditions to us and vs
        call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        
        end do

        ! Form the RHS of the pressure poisson equation
        call calculate_rhs(M,us,vs,R,rho,dt)

        ! Solve for pressure
        call calculate_pressure_amgx(A,P,R,init_status)

        ! Perform the corrector step to obtain the velocity
        call corrector(M,u,v,us,vs,P,rho,dt)

        ! Initialize velocity to zero and then Interpolate velocity
        call cilia(:)%set_U(0.0d0)
        call particles%set_velocity(0.0d0)
        call interpolate_velocity(M,cilia,u,v)
        call interpolate_velocity(M,particles,u,v)

        call vorticity(M,u,v,w)
        call cilia(:)%set_mden(0.0d0)
        call interpolate_vorticity(M,cilia,w)

        ! Update structure
        call cilia(:)%update(dt,1,1)
        call particles(:)%update_position(dt)

        ! update "guess value" of pressure
        ! Pold=p+Pold

        ! Write files every Nth timestep
        if (mod(it,it_save).eq.0) then 
        ! if (mod(it,1).eq.0) then 
            call write_field(u,'u',it) 
            call write_field(v,'v',it)
            call write_field(w,'w',it); 
            do ic = 1,ncilia 
                call write_field(cilia(ic)%XE,'C',it)
            end do
            call write_field(particles(1)%XE,'P',it)
        end if
    end do

end program ibmc
