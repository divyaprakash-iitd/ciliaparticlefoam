module soft_cilia
    ! List of Magic Numbers
    ! Domain size (Lx,Ly,Lz)
    ! Z-coordinate to convert 2d case to 3d

    ! use fem3d
    use mod_cilia
    use mod_io
    implicit none

    ! Parameters
    real(8), parameter :: PI = 3.141592653589793
    
    ! Computational Domain
    real(8)    :: Lx = 0.01 !50 !0.3
    real(8)    :: Ly = 0.002 !50 !0.1
    real(8)    :: Lz = 2.5e-5 !50 !0.1
    
    ! ! Particle information
    ! integer         :: nvp ! Number of vertices in the particle
    ! real(real64)    :: aa, bb ! Major and minor axis
    ! real(8)         :: Kp, Bp
    ! integer(int32)  :: itnum ! iteration number
    
    ! Cilia information
    real(8) :: x0(2)
    real(8) :: ftip  
    real(8) :: mtip  
    integer :: ncilia
    integer :: icilia
    real(8) :: l
    integer :: nvc ! Number of vertices in the cilia
    type(cilium), allocatable :: cilia(:)
    real(8) :: Kc, Bc, kap, bap
    real(8), allocatable :: xcil(:,:)
    real(8) :: lcbed, delc
    logical :: orient, ctype
    real(8) :: nplus, nminus

    ! Pseudo 3D
    real(8) :: zcoord !Since cilia is ony 2D, we assign it's z-coordinate manually

    ! Keeping track of iterations for writing purposes
    integer(int32)  :: itnum ! iteration number

    ! FEM data
    ! type(festruct), allocatable :: particles(:)

    ! Namelists for input
    ! namelist /particleprops/ Kp, nvp, Bp, aa, bb
    namelist /ciliaprops/ l, ftip, mtip, ncilia, Kc, nvc, Bc, kap, bap

    !---------------------- Begin Calculations ------------------------------------!

    contains 

    subroutine say_hello() bind(C)
       use iso_c_binding, only: C_INT, C_CHAR
       implicit none
       print *, "Hello!"
    end subroutine say_hello

    subroutine generatecilia(noelpts,cdl,h) bind(C)
    use iso_c_binding, only: C_INT, C_CHAR, c_double
    implicit none

    integer(C_INT), intent(inout) :: noelpts ! No. of points in all the cilia
    real(c_double), intent(in) :: h ! Mesh width in openfoam
    real(c_double), intent(inout) :: cdl ! The segment length of cilia

    integer(int32) :: i, err

    ! Read input data from file
    open(1004,file="input_params.dat",form='formatted')
    READ(unit=1004,nml=ciliaprops,iostat=err)
    close(1004)
    
    ! Copy the total no. of points to the openfoam variable
    noelpts = ncilia*(nvc-1) ! (-1) because the cilia is only represented by midpoints
    print *, "ncilia = ", noelpts
    
    ! Create cilia
    x0 = [Lx/4.0, 3*h]
    ! l = Ly/4.0d0 
    allocate(cilia(ncilia), xcil(ncilia,2))
    lcbed = 0.80d0 * Lx ! Length of cilia bed
    delc = lcbed / (ncilia - 1)
    xcil (:,2) = x0(2)
    xcil(:,1) = (Lx - lcbed)/2.0d0 + delc * [(i-1, i=1, ncilia)]
    ! xcil(:,1) = Lx/2.0

    ! Detection cilia
    orient = .FALSE.
    ctype = .FALSE.
    nplus = 0.0d0
    nminus = 0.0d0
    do i = 1,ncilia
        cilia(i) = cilium(xcil(i,:),l,nvc,Kc,Bc,kap,bap,orient,ctype,nplus,nminus)
    end do
            
    do i = 1,ncilia 
        call write_field(cilia(i)%XE,'C',1)
    end do

    ! Initialize the iteration number
    itnum = 1
    cdl = cilia(1)%dl

    end subroutine generatecilia
   
    ! subroutine arraycheck(pxyz,n) bind(C)
    !    use iso_c_binding, only: c_int, c_double, c_loc
    !    implicit none
       
    !    integer(c_int), intent(in) :: n
    !    real(c_double), intent(inout)   :: pxyz(n)
    !    ! integer(c_int), intent(inout)   :: pxyz(:)
    !    ! integer(c_int), intent(inout)   :: pxyz(n)
       
    !    integer(int32) :: i, npoints
    !    npoints = size(particles(1)%XE,1)
      
    !    ! print *, "SIZE: ", size(pxyz)
    !    ! print *, "pxyz1: ", pxyz(5)
    !    do i = 1,5!npoints
    !        ! pxyz(i)   = particles(1)%XE(i,1)
    !        pxyz(i)   = 23651491.23
    !    end do


    ! end subroutine arraycheck

    subroutine getpositions(XC,YC,ZC,nn) bind(C)
       ! It takes in the position arrays defined in openfoam and fills
       ! it with the particle's position values
       use iso_c_binding, only: c_int, c_double, c_loc
       implicit none

       integer(c_int), intent(in)      :: nn
       real(c_double), intent(inout)   :: XC(nn),YC(nn),ZC(nn)

       ! Create a matrix to store all the cilia nodes
       real(c_double), allocatable :: carrayx(:,:), carrayy(:,:), carrayz(:,:)
       integer(int32) :: i
       
       allocate(carrayx(nvc-1,ncilia),carrayy(nvc-1,ncilia),carrayz(nvc-1,ncilia))

       !Since cilia is ony 2D, we assign it's z-coordinate manually
       zcoord = Lz/2.0d0
       do i = 1,ncilia
           carrayx(:,i) = cilia(i)%XE(:,1)
           carrayy(:,i) = cilia(i)%XE(:,2)
           carrayz(:,i) = zcoord !cilia(i)%XE(:,3)
       end do
    
       ! Reshape and copy to openfoam array
        XC   = reshape(carrayx,[nn])
        YC   = reshape(carrayy,[nn])
        ZC   = reshape(carrayz,[nn])
       
    end subroutine getpositions
   
    subroutine calculateforcesandmoments(FXC,FYC,FZC,MXC,MYC,MZC,nn) bind(C)
       ! Calculates the forces in the particle
       ! Transfers those forces to the arrays passed in by openfoam
       ! Addition: Calculate moments in the cilia and pass those to the openfoam arrays as well
       use iso_c_binding, only: c_int, c_double, c_loc
       implicit none

       integer(c_int), intent(in)      :: nn
       real(c_double), intent(inout)   :: FXC(nn),FYC(nn),FZC(nn)
       real(c_double), intent(inout)   :: MXC(nn),MYC(nn),MZC(nn)

       integer(int32) :: i

       ! Create a matrix to store all the cilia forces and moments
       real(c_double), allocatable :: cfx(:,:), cfy(:,:), cfz(:,:)
       real(c_double), allocatable :: cmx(:,:), cmy(:,:), cmz(:,:)
       
       allocate(cfx(nvc-1,ncilia),cfy(nvc-1,ncilia),cfz(nvc-1,ncilia))
       allocate(cmx(nvc-1,ncilia),cmy(nvc-1,ncilia),cmz(nvc-1,ncilia))

       ! Add lines here to calculate the moments as well which will be copied to the
       ! array in openfoam
       do i = 1,ncilia
           call cilia(i)%forces(0.000000000001d0)
           call cilia(i)%moments()
       end do

       ! Gather forces
       do i = 1,ncilia
           cfx(:,i) = cilia(i)%fden(:,1)
           cfy(:,i) = cilia(i)%fden(:,2)
           cfz(:,i) = 0.0d0 !cilia(i)%fden(:,3)
       end do
       ! Copy to openfoam
       FXC = reshape(cfx,[nn])
       FYC = reshape(cfy,[nn])
       FZC = reshape(cfz,[nn])
       
       ! Gather moments
       do i = 1,ncilia
           cmx(:,i) = 0.0d0 !cilia(i)%mden(:,1)
           cmy(:,i) = 0.0d0 !cilia(i)%mden(:,2)
           cmz(:,i) = cilia(i)%mden !cilia(i)%mden(:,1)
       end do
       ! Copy to openfoam
       MXC = reshape(cmx,[nn])
       MYC = reshape(cmy,[nn])
       MZC = reshape(cmz,[nn])

    end subroutine calculateforcesandmoments


    subroutine updatepositions(U,V,W,MX,MY,MZ,dt,nn) bind(C)
       ! Take in the velocity and angular velocity from openfoam
       ! Use them to update the cilia position and orientation
       use iso_c_binding, only: c_int, c_double, c_loc
       implicit none

       integer(c_int), intent(in)      :: nn
       real(c_double), intent(inout)   :: U(nn), V(nn) ,W(nn)
       real(c_double), intent(inout)   :: MX(nn), MY(nn) ,MZ(nn)
       real(c_double), intent(in)      :: dt

       integer(int32) :: i

       ! Create a matrix to store all the interpolated cilia velocity
       real(c_double), allocatable :: cux(:,:), cuy(:,:), cuz(:,:)
       ! Create a matrix to store all the interpolated cilia moments
       real(c_double), allocatable :: cmx(:,:), cmy(:,:), cmz(:,:)
       
       allocate(cux(nvc-1,ncilia),cuy(nvc-1,ncilia),cuz(nvc-1,ncilia))
       allocate(cmx(nvc-1,ncilia),cmy(nvc-1,ncilia),cmz(nvc-1,ncilia))
       
       ! Rehshape the obtained velocities in a matrix with columns corresponding to each cilia
        cux = reshape(U,[nvc-1,ncilia])
        cuy = reshape(V,[nvc-1,ncilia])
        cuz = reshape(W,[nvc-1,ncilia])
       
        ! Rehshape the obtained angular velocities in a matrix with columns corresponding to each cilia
        cmx = reshape(MX,[nvc-1,ncilia])
        cmy = reshape(MY,[nvc-1,ncilia])
        cmz = reshape(MZ,[nvc-1,ncilia])
   
        ! Copy openfoam data to cilia 
       do i = 1,ncilia
           cilia(i)%U(:,1) = cux(:,i)
           cilia(i)%U(:,2) = cuy(:,i)
        !    cilia(i)%U(:,3) = cuz(:,i)
       end do
       
       do i = 1,ncilia
        !    cilia(i)%mden(:,1) = cmx(:,i)
        !    cilia(i)%mden(:,2) = cmy(:,i)
           cilia(i)%mden = cmz(:,i)
       end do

       itnum = itnum + 1
       
       if (mod(itnum,200).eq.0) then
            do i = 1,ncilia
                call write_field(cilia(i)%XE,'C',itnum)
            end do
       end if

       do i = 1,ncilia
           call cilia(i)%update(dt)
       end do


        do i = 1,ncilia 
            call write_field(cilia(i)%XE,'L',1)
        end do
    end subroutine updatepositions

    ! ! Create 3 subroutines
    ! ! 1. Creates the ellipse and it's coordinates and connectivity. Basically reads it form the python generated file.
    ! ! 2. Calls the force calculation and fills up the force vector.
    ! ! 3. Takes the velocity vector from the C program and uses it to update the particle's nodes positions.

end module soft_cilia
