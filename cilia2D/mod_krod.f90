module mod_krod
    use iso_fortran_env, only: int32, real64
    use mod_solid
    implicit none

    ! real(real64), parameter :: PI = 3.141592653589793
    real(real64), parameter :: PI=4.D0*DATAN(1.0D0)

    private :: PI

    type, extends(solid) :: krod
    ! type, abstract :: krod
        integer(int32) :: NV ! Number of vertices
        integer(int32) :: NE ! Number of elements
        real(real64):: X0(2) ! Origin of the cilium
        real(real64) :: L ! Length of cilium
        real(real64) :: K ! same as b1 and b3 
        real(real64) :: B ! same as a2 
        ! real(real64) :: dl 
        real(real64), allocatable :: XV(:,:) ! Coordinates of the vertices
        real(real64), allocatable :: XEOLD(:,:) ! Coordinates of the mid-points
        ! real(real64), allocatable :: XE(:,:) ! Coordinates of the mid-points
        real(real64), allocatable :: F(:,:) ! Internal force defined at vertices
        ! real(real64), allocatable :: U(:,:) ! Velocity defined at midpoint
        real(real64), allocatable :: Uold(:,:) ! Velocity defined at midpoint
        real(real64), allocatable :: Umid(:,:) ! Velocity defined at midpoint
        real(real64), allocatable :: that(:,:) ! tangent vector defined at midpoint
        real(real64), allocatable :: nhat(:,:) ! normal vector defined at midpoint
        real(real64), allocatable :: thatv(:,:) ! tangent vector defined at vertex 
        real(real64), allocatable :: nhatv(:,:) ! normal vector defined at vertex
        real(real64), allocatable :: sigma(:) ! Tangential force defined at vertices
        real(real64), allocatable :: N(:) ! Normal force defined at vertices
        real(real64), allocatable :: Mom(:) ! Internal moment defined at vertices
        ! real(real64), allocatable :: fden(:,:) ! force density on midpoint
        real(real64), allocatable :: mden(:) ! moment density defined at midpoint
        real(real64), allocatable :: mdenold(:) ! moment density defined at midpoint
        real(real64), allocatable :: mdenmid(:) ! moment density defined at midpoint
        real(real64), allocatable :: phi(:) ! Angle at element midpoints
        real(real64), allocatable :: phiold(:) ! Angle at element midpoints
        real(real64), allocatable :: phiv(:) ! Angle at vertices
        
        real(real64), allocatable :: ds(:) ! element length
        ! real(real64) :: rada ! Actual radius
        real(real64) :: kap ! Anchor point linear spring constant
        real(real64) :: bap ! Anchor point angular spring constant
        real(real64) :: XM0(2) ! Original anchor location
        real(real64), allocatable :: K0(:) ! Reference curvature
        real(real64) :: phi0org ! Angle required for moment anchor
        real(real64), allocatable :: thetav(:) ! Ellipse angles at the vertices
        real(real64), allocatable :: theta(:) ! Ellipse angles at the mid-points
        real(real64), allocatable :: sigma0(:), N0(:) ! Reference stretching forces

        ! Add variables to sliding displacment, sliding velocity and motor force
        ! Remember these are scalar as they are defined along the cilia
        real(real64), allocatable :: SD(:) ! Sliding displacement defined at midpoint
        real(real64), allocatable :: SV(:) ! Sliding velocity defined at midpoint
        real(real64), allocatable :: motorf(:) ! motor force defined at midpoint
        real(real64), allocatable :: nplus(:)  ! fraction of motors in the bound state on the "+" filament per segment
        real(real64), allocatable :: nminus(:) ! fraction of motors in the bound state on the "-" filament per segment
        logical :: orient, release, ctype 

        
    end type krod 

end module mod_krod
