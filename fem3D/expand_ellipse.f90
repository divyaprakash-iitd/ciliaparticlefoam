program main
    use fem3d
    use iso_fortran_env, only: int32, real64
    implicit none

    character(len=256) :: file_path
    real(real64), allocatable :: M(:,:), XE(:,:), FN(:,:), UN(:,:), rmag(:), XORG(:,:)
    type(festruct) :: particle
    integer(int32) :: i, j, npoints, niter, iter
    real(real64) :: fmag, dt, mu
    real(real64) :: co, kval, dl


    file_path = "data/M_file.bin"
    call read_data(file_path,M)
    
    file_path = "data/XE_file.bin"
    call read_data(file_path,XE)

    FN = 0.0d0*XE
    UN = 0.0d0*XE

    ! Calculate the shapecoefficients
    co = 5000.0d0
    kval = 10000.0d0
    dl = 1.0d0
    particle = festruct(int(M),XE,FN,UN,co,kval,dl) ! kp = kval, co = bp , dl = 1.0d0

    ! Apply radial forces
    ! Calculate the distance from the center of each node
    npoints = size(XE,1)
    allocate(rmag(npoints),XORG(npoints,3))
    do i = 1,npoints
        rmag(i) = norm2(XE(i,:))
    end do
    XORG = XE


    ! print *, shape(XORG)
    niter = 500
    dt = 0.001
    mu = 1000
    fmag = 5000

    do iter = 1,niter 
        ! Initialize the forces at the ellipse points to be zero
        particle%fden = 0.0d0
        ! Shape of fden (npoints,3)
        ! Apply the expansion forces
        if (iter.lt.50) then
            mu = 1000
            do i = 1,npoints
                if (rmag(i).ge.0.95d0) then
                    particle%fden(i,:) = fmag * XORG(i,:)/rmag(i)
                    ! print *, "fmag=", fmag
                end if
            end do
        else 
            mu = 100
        end if
        
        call particle%calculate_forces()

        ! Update the point's locations
        do i = 1,npoints
            particle%XE(i,:) = particle%XE(i,:) + dt/mu * particle%fden(i,:)
        end do

        if (mod(iter,5).eq.0) then
            print *, "iter = ", iter
            call write_field(particle%XE,'X',iter)
            call write_field(particle%fden,'F',iter)
        end if
    end do

end program main
