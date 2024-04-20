module soft_cilia
    ! use fem3d
    use mod_cilia
    use mod_io
    implicit none

    ! Parameters
    real(8), parameter :: PI = 3.141592653589793
    
    ! Computational Domain
    real(8)    :: Lx = 10 !50 !0.3
    real(8)    :: Ly = 10 !50 !0.1
    real(8)    :: Lz = 10 !50 !0.1
    
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

    subroutine generatecilia(noelpts,h) bind(C)
    use iso_c_binding, only: C_INT, C_CHAR, c_double
    implicit none

    integer(C_INT), intent(inout) :: noelpts ! No. of points in all the cilia
    real(c_double), intent(in) :: h ! Mesh width in openfoam

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
    l = Ly/4.0d0 
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

    ! subroutine getpositions(XC,YC,ZC,nn) bind(C)
    !    ! It takes in the position arrays defined in openfoam and fills
    !    ! it with the particle's position values
    !    use iso_c_binding, only: c_int, c_double, c_loc
    !    implicit none

    !    integer(c_int), intent(in)      :: nn
    !    real(c_double), intent(inout)   :: XC(nn),YC(nn),ZC(nn)

    !    integer(int32) :: i, nparticles, npoints

    !    ! print *, "Size of XC: ", size(XC)
    !    nparticles = 1

    !    npoints = size(particles(1)%XE,1)

    !    do i = 1,npoints
    !        XC(i)   = particles(1)%XE(i,1)
    !        YC(i)   = particles(1)%XE(i,2)
    !        ZC(i)   = particles(1)%XE(i,3)
    !    end do
    ! end subroutine getpositions
   
    ! subroutine calculateforces(FXC,FYC,FZC,nn) bind(C)
    !    ! Calculates the forces in the particle
    !    ! Transfers those forces to the arrays passed in by openfoam
    !    ! Addition: Calculate moments in the cilia and pass those to the openfoam arrays as well
    !    use iso_c_binding, only: c_int, c_double, c_loc
    !    implicit none

    !    integer(c_int), intent(in)      :: nn
    !    real(c_double), intent(inout)   :: FXC(nn),FYC(nn),FZC(nn)

    !    integer(int32) :: i, nparticles, npoints

    !    nparticles = 1

    !    npoints = size(particles(1)%XE,1)


    !    ! Add lines here to calculate the moments as well which will be copied to the
    !    ! array in openfoam
    !    do i = 1,nparticles
    !        call particles(i)%calculate_forces()
    !    end do

    !    do i = 1,npoints
    !        FXC(i)  = particles(1)%fden(i,1)
    !        FYC(i)  = particles(1)%fden(i,2)
    !        FZC(i)  = particles(1)%fden(i,3)
    !    end do
    ! end subroutine calculateforces


    ! subroutine updatepositions(U,V,W,dt,nn) bind(C)
    !    ! Take in the velocity and angular velocity from openfoam
    !    ! Use them to update the cilia position and orientation
    !    use iso_c_binding, only: c_int, c_double, c_loc
    !    implicit none

    !    integer(c_int), intent(in)      :: nn
    !    real(c_double), intent(inout)   :: U(nn), V(nn) ,W(nn)
    !    real(c_double), intent(in)      :: dt

    !    integer(int32) :: i, nparticles, npoints

    !    nparticles = 1

    !    npoints = size(particles(1)%XE,1)
       
    !    do i = 1,npoints
    !        particles(1)%U(i,1) = U(i)
    !        particles(1)%U(i,2) = V(i)
    !        particles(1)%U(i,3) = W(i)
    !    end do

    !    itnum = itnum + 1
       
    !    if (mod(itnum,200).eq.0) then
    !        call write_field(particles(1)%XE,'P',itnum)
    !    end if

    !    do i = 1,nparticles
    !        call particles(i)%update_position(dt)
    !    end do

    ! end subroutine updatepositions

    ! ! Create 3 subroutines
    ! ! 1. Creates the ellipse and it's coordinates and connectivity. Basically reads it form the python generated file.
    ! ! 2. Calls the force calculation and fills up the force vector.
    ! ! 3. Takes the velocity vector from the C program and uses it to update the particle's nodes positions.

end module soft_cilia
