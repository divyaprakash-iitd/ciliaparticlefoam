module mod_cilia
    use iso_fortran_env, only: int32, real64
    use mod_krod
    implicit none
    
    ! real(real64), parameter :: PI = 3.141592653589793
    real(real64), parameter :: PI=4.D0*DATAN(1.0D0)

    private :: PI
    
    type, extends(krod) :: cilium
        ! integer(int32) :: NV ! Number of vertices
        ! integer(int32) :: NE ! Number of elements
        ! real(real64):: X0(2) ! Origin of the cilium
        ! real(real64) :: L ! Length of cilium
        ! real(real64) :: K ! same as b1 and b3 
        ! real(real64) :: B ! same as a2 
        ! real(real64) :: dl 
        ! real(real64), allocatable :: XV(:,:) ! Coordinates of the vertices
        ! real(real64), allocatable :: XEOLD(:,:) ! Coordinates of the mid-points
        ! real(real64), allocatable :: XE(:,:) ! Coordinates of the mid-points
        ! real(real64), allocatable :: F(:,:) ! Internal force defined at vertices
        ! real(real64), allocatable :: U(:,:) ! Velocity defined at midpoint
        ! real(real64), allocatable :: that(:,:) ! tangent vector defined at midpoint
        ! real(real64), allocatable :: nhat(:,:) ! normal vector defined at midpoint
        ! real(real64), allocatable :: thatv(:,:) ! tangent vector defined at vertex 
        ! real(real64), allocatable :: nhatv(:,:) ! normal vector defined at vertex
        ! real(real64), allocatable :: sigma(:) ! Tangential force defined at vertices
        ! real(real64), allocatable :: N(:) ! Normal force defined at vertices
        ! real(real64), allocatable :: Mom(:) ! Internal moment defined at vertices
        ! real(real64), allocatable :: fden(:,:) ! force density on midpoint
        ! real(real64), allocatable :: mden(:) ! moment density defined at midpoint
        ! real(real64), allocatable :: phi(:) ! Angle at element midpoints
        ! real(real64), allocatable :: phiold(:) ! Angle at element midpoints
        ! real(real64), allocatable :: phiv(:) ! Angle at vertices
        
        ! real(real64), allocatable :: ds(:) ! element length
        ! ! real(real64) :: rada ! Actual radius
        ! real(real64) :: kap ! Anchor point linear spring constant
        ! real(real64) :: bap ! Anchor point angular spring constant
        ! real(real64) :: XM0(2) ! Original anchor location
        
        contains
            procedure :: unit_vectors
            procedure :: forces
            procedure :: moments
            procedure :: velocity
            procedure :: update
            procedure :: set_U
            procedure :: set_mden

    end type cilium 

    interface cilium
        module procedure :: cilium_constructor
    end interface cilium

    contains 
    ! Write a constructor which only requires number of points and fills in the
    ! other information by itself

    pure type(cilium) function cilium_constructor(X0, L, NV, K, B, kap, bap) result(self)
       real(real64), intent(in) :: X0(2) ! Length of cilium
       real(real64), intent(in) :: L ! Length of cilium
       integer(int32), intent(in) :: NV       
       real(real64), intent(in) :: K
       real(real64), intent(in) :: B
       real(real64), intent(in) :: kap
       real(real64), intent(in) :: bap

       integer(int32) :: i
       real(real64) :: theta, theta1, theta2, dtheta
       real(real64) :: radius

       self%NV = NV
       self%L = L
       self%X0 = X0

       self%NE = NV - 1
       self%K = K
       self%B = B
       self%kap = kap
       self%bap = bap
       
       ! Allocate the arrays
       ! Defined at vertices
       allocate(self%XV(self%NV,2), self%thatv(self%NV,2), &
                    self%nhatv(self%NV,2), self%phiv(self%NV), &
                    self%sigma(self%NV), self%N(self%NV), &
                    self%F(self%NV,2), self%Mom(self%NV))
       
       ! Defined at mid-points
       allocate(self%U(self%NE,2), self%that(self%NE,2), self%nhat(self%NE,2), &
                     self%XE(self%NE,2), self%phi(self%NE), & 
                     self%ds(self%NE-1), self%fden(self%NE,2), self%mden(self%NE), &
                     self%xeold(self%NE,2), self%phiold(self%NE), &
                     self%Uold(self%NE,2), self%mdenold(self%NE), &
                     self%Umid(self%NE,2), self%mdenmid(self%NE))

       ! Initialize the variables
       self%XV = 0.0d0
       self%that = 0.0d0
       self%nhat = 0.0d0
       self%thatv = 0.0d0
       self%nhatv = 0.0d0
       self%phi = 0.0d0!PI/2.0d0
       self%phiv = 0.0d0!PI/2.0d0
       self%sigma = 0.0d0
       self%N = 0.0d0
       self%F = 0.0d0
       self%U = 0.0d0
       self%XE = 0.0d0
       self%ds = 0.0d0
       self%Mom = 0.0d0
       self%fden =0.0d0
       self%mden = 0.0d0
       self%Uold = 0.0d0
       self%mdenold = 0.0d0
       self%Umid = 0.0d0
       self%mdenmid = 0.0d0

    ! !-------------------------Horizontal Cilia-----------------------------!
    !    ! Calculate the vertices of the cilium
    !    self%dl = L/self%NE
    !    self%XV(:,2) = X0(2)
    !    self%XV(:,1) = X0(1) + self%dl * [(i, i=0, self%NE)]
      
    !    ! Calculate the mid-points of the elements 
    !    self%XE = self%XV(1:self%NE,:) + (self%XV(2:self%NV,:) - self%XV(1:self%NE,:)) / 2.0d0
       
    ! !------Add code to calculate unit vectors for the first time step------!
    !    ! In the first time step, all the cilia are horizontal
        ! self%phi = 0.0d0
        ! self%phiv = 0.0d0
    ! !-----------------------------------------------------------------------!

    !-------------------------Vertical Cilia-----------------------------!
       ! Calculate the vertices of the cilium
       self%dl = L/self%NE
       self%XV(:,1) = X0(1)
       self%XV(:,2) = X0(2) + self%dl * [(i, i=0, self%NE)]
      
       ! Calculate the mid-points of the elements 
       self%XE = self%XV(1:self%NE,:) + (self%XV(2:self%NV,:) - self%XV(1:self%NE,:)) / 2.0d0
       self%XEOLD = self%XE
        ! Original anchor point location
        self%XM0 = self%XE(1,:)
    !------Add code to calculate unit vectors for the first time step------!
       ! In the first time step, all the cilia are vertical
        self%phi = PI/2.0d0
        self%phiv = PI/2.0d0
    !-----------------------------------------------------------------------!

    ! !------------------------- Circular Cilia ------------------------!       
    !    ! Calculate the mid-points of the cilium
    !    radius = 2*L
    !    theta1 = PI/3.0d0
    !    theta2 = 2.0d0/3.0d0 * PI
    !    dtheta = (theta2 - theta1)/self%NE
    !    theta = theta1 + dtheta/2.0d0
    !    do i = 1,self%NE
    !         self%XE(i,1) = radius * cos(theta) 
    !         self%XE(i,2) = radius * sin(theta) 
    !         theta = theta + dtheta
    !    end do

    !    ! Calculate the vertices
    !    do i = 2,self%NV-1
    !         self%XV(i,1) = 0.5d0 * (self%XE(i-1,1) + self%XE(i,1))
    !         self%XV(i,2) = 0.5d0 * (self%XE(i-1,2) + self%XE(i,2))
    !    end do

    !    ! Calculate the first and last vertex
    !    self%XV(1,:) = 0.5d0 * (self%XE(1,:) - [radius*cos(theta1-dtheta/2.0d0), radius*sin(theta1-dtheta/2.0d0)])
    !    self%XV(self%NV,:) = 0.5d0 * (self%XE(self%NE,:) - [radius*cos(theta2+dtheta/2.0d0), radius*sin(theta2+dtheta/2.0d0)])

    !    ! Translate the points to cilia origin
    !    self%XE(:,1) = self%XE(:,1) + self%X0(1) 
    !    self%XE(:,2) = self%XE(:,2) + self%X0(2) - radius
    !    self%XV(:,1) = self%XV(:,1) + self%X0(1) - radius
    !    self%XV(:,2) = self%XV(:,2) + self%X0(2) 

    !    ! Element length
    ! !    self%dl = norm2([self%XV(2,1) - self%XV(1,1), self%XV(2,2)- self%XV(1,2)])
    !     self%dl = sqrt((self%XE(3,1)-self%XE(2,1))**2 + (self%XE(3,2)-self%XE(2,2))**2)

    !    ! Initialize the phi and phiv
    !    ! Phi is the angle that the tangent makes with the x-axis defined at mid-points
    !    self%phi(1) = theta1 + dtheta/2.0d0
    !    do i = 2,self%NE
    !         self%phi(i) = self%phi(i-1) + dtheta 
    !    end do
    !    self%phi = self%phi + PI/2.0d0 
    !    self%phiv(1) = (theta1-dtheta/2.0d0) + PI/2.0d0
    !    self%phiv(self%NV) = (theta2+dtheta/2.0d0) + PI/2.0d0 ! The rest are calculated in the unit vectors subroutine

    !    self%rada = norm2(self%XE(1,:)-self%X0)
    ! !--------------------------------------------------------------------!
        
    end function cilium_constructor

    pure subroutine unit_vectors(self)
        class(cilium), intent(inout) :: self

        self%phiv(2:self%NV-1)=(0.5d0)*(self%phi(1:self%NE-1) + self%phi(2:self%NE))

        ! Calculate the tangential vectors
        self%that(1:self%NE,1) = cos(self%phi)
        self%that(1:self%NE,2) = sin(self%phi)
        
        self%nhat(1:self%NE,1) = -sin(self%phi)
        self%nhat(1:self%NE,2) = cos(self%phi)

        self%thatv(2:self%NV-1,1) = cos(self%phiv(2:self%NV-1))
        self%thatv(2:self%NV-1,2) = sin(self%phiv(2:self%NV-1))
        
        self%nhatv(2:self%NV-1,1) = -sin(self%phiv(2:self%NV-1))
        self%nhatv(2:self%NV-1,2) =  cos(self%phiv(2:self%NV-1))
     
    end subroutine unit_vectors

    elemental impure subroutine moments(self,mtip)
        class(cilium), intent(inout) :: self
        real(real64), optional, intent(in) :: mtip
        
        integer(int32) :: k
        real(real64) :: part1(self%NE)
        real(real64) :: part2(self%NE)
        real(real64) :: momap, theta

        part1 = 0.0d0 ! Initialize
        part2 = 0.0d0 ! Initialize

        ! Initialize mom and mden everytime it's called
        self%mom = 0.0d0
        self%mden = 0.0d0
 
        
        !------Moments at the vertices------!
        ! Let us take ds = dl
        ! The moments at the first and last vertex is zero
        do k = 2,self%NV-1
            self%mom(k) = (self%B / self%dl) * DOT_PRODUCT((self%that(k,:)-self%that(k-1,:)), & 
                                                     self%nhatv(k,:))
        end do

        !------Moment density at the mid-points------!
        self%mden = ( self%mom(2:self%NV) - self%mom(1:self%NV-1) ) / self%dl

        ! print *, self%mom
        ! print *, "------------------------------"
        ! print *, ( self%mom(2:self%NV) - self%mom(1:self%NV-1) )

        do k = 1,self%NE

            if (k.eq.self%NE) then
                part1(k) = 0.0d0
            else
                part1(k) = (self%XE(k+1,1) - self%XE(k,1)) * self%F(k+1,2) - &
                           (self%XE(k+1,2) - self%XE(k,2)) * self%F(k+1,1)
            end if

            if (k.eq.1) then
                part2(k) = 0.0d0
            else
                part2(k) = (self%XE(k,1) - self%XE(k-1,1)) * self%F(k,2) - &
                           (self%XE(k,2) - self%XE(k-1,2)) * self%F(k,1)
            end if
        end do

        self%mden = self%mden + 0.5d0 * (part1 + part2) / self%dl

        ! Apply anchor point moment
        theta = PI/2.0d0 - atan2(self%that(1,2), self%that(1,1))
        !print *, "THe = ", theta
        
        momap = -self%bap * theta
        self%mden(1) = self%mden(1) + (-momap)
        
        ! Add moments to the tips
        if (present(mtip)) then
            self%mden(self%NE) = self%mden(self%NE) + mtip
        end if

        ! print *, sum(part1)
        ! print *, sum(part2)
        !print *, momap, theta, self%that(1,2), self%that(1,1)
        

    end subroutine moments

    elemental impure subroutine forces(self,ftip)
        class(cilium), intent(inout) :: self
        real(real64), optional, intent(in) :: ftip

        INTEGER :: i
        real(real64) :: dap ! Distance from anchor point
        real(real64) :: fap(2) ! Anchorage force
        real(real64), allocatable :: dxds(:,:)
        
        allocate(dxds(1:self%NV,2))
        ! Calculate the unit vectors
        call self%unit_vectors()
                
        dxds=0.0d0;        
        dxds(2:self%NV-1,:) = (self%XE(2:self%NE,:) - self%XE(1:self%NE-1,:))/self%dl

        ! do i = 1,self%NV 
        !     print *, DOT_PRODUCT(dxds(i,:),self%nhatv(i,:))
        !     print *, DOT_PRODUCT(dxds(i,:),self%thatv(i,:))
        ! end do

        ! Initialize forces everytime it's called
        self%sigma = 0.0d0
        self%N = 0.0d0
        self%F = 0.0d0
        self%fden = 0.0d0

        ! Calculate tangential force
        self%sigma(2:self%NV-1) = self%K * (dxds(2:self%NV-1,1)*self%thatv(2:self%NV-1,1) + &
					   dxds(2:self%NV-1,2)*self%thatv(2:self%NV-1,2) - 1.0d0)

        ! Calculate normal force
        self%N(2:self%NV-1) = self%K * (dxds(2:self%NV-1,1)*self%nhatv(2:self%NV-1,1) + &
					   dxds(2:self%NV-1,2)*self%nhatv(2:self%NV-1,2))
 
       ! Calculate the forces
       ! AB: Make the normal force zero
       self%F(:,1) = self%sigma * self%thatv(:,1) + self%N * self%nhatv(:,1)
       self%F(:,2) = self%sigma * self%thatv(:,2) + self%N * self%nhatv(:,2)
        

        ! if (present(ftip)) then
            ! self%F(self%NV,1)  = self%F(self%NV,1) + ftip
        ! end if

        ! print *, self%sigma 
        ! print *, self%N 
       deallocate(dxds)
        
       
       !------Force density at the mid-points------!
       self%fden = (1.0d0 / self%dl) * ( self%F(2:self%NV,:) - self%F(1:self%NV-1,:) )

        if (present(ftip)) then
            self%fden(self%NE,1)  = self%fden(self%NE,1) + ftip
            ! self%fden(self%NE,1)  = self%fden(1,2) - ftip
            ! self%fden(self%NE,2) = self%fden(self%NE,2) + ftip
            ! self%F(1,2)  = self%F(1,2) + ftip
            ! self%F(self%NV,2) = self%F(self%NV,2) + ftip
        end if
        
        ! Apply anchor point force
        ! Calculate the distance of the first mid-point from the anchor point
        dap = norm2([(self%XE(1,1)-self%X0(1)), (self%XE(1,2)-self%X0(2))])
       
        fap(1) = self%kap * (self%XM0(1)-self%XE(1,1))
        fap(2) = self%kap * (self%XM0(2)-self%XE(1,2))
        !print *, "Fap = ", fap
        
        self%fden(1,:)  = self%fden(1,:) + fap

        ! print *, "maxfden1 = ", maxval(self%fden(:,1))
        ! print *, "maxfden2 = ", maxval(self%fden(:,2))

    end subroutine forces

    subroutine velocity(self)
        class(cilium), intent(inout) :: self

        self%U = ((self%F(2:self%NV,:) - self%F(1:self%NV-1,:))) / self%dl
    end subroutine velocity
    
    elemental subroutine update(self, dt, den, it)
        class(cilium), intent(inout) :: self
        integer(int32), intent(in) :: den
        real(real64), intent(in) :: dt
        integer(int32), intent(in) :: it
        integer j

        self%Umid(1:self%NE,1) = 0.50d0 * (self%U(1:self%NE,1) + self%Uold(1:self%NE,1))
        self%Umid(1:self%NE,2) = 0.50d0 * (self%U(1:self%NE,2) + self%Uold(1:self%NE,2))
        self%mdenmid(1:self%NE) = 0.50d0 * (self%mden(1:self%NE) + self%mdenold(1:self%NE))

        ! if (it.eq.1) then
            self%Umid = self%U
            self%mdenmid = self%mden
        ! end if
        
        self%XE(1:self%NE,1) = self%XE(1:self%NE,1) + dt * self%Umid(1:self%NE,1)
        self%XE(1:self%NE,2) = self%XE(1:self%NE,2) + dt * self%Umid(1:self%NE,2)

        !------Update phi at the mid-points------!
        self%phi(1:self%NE) = self%phi(1:self%NE) + dt * self%mdenmid(1:self%NE)
    end subroutine update

    elemental subroutine set_U(self, val)
        class(cilium), intent(inout) :: self
        real(real64), intent(in) :: val

        self%U = val
    end subroutine set_U
    
    elemental subroutine set_mden(self, val)
        class(cilium), intent(inout) :: self
        real(real64), intent(in) :: val

        self%mden = val
    end subroutine set_mden

end module mod_cilia    
