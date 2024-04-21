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
            procedure :: calculate_motorf

    end type cilium 

    interface cilium
        module procedure :: cilium_constructor
    end interface cilium

    contains 
    ! Write a constructor which only requires number of points and fills in the
    ! other information by itself

    pure type(cilium) function cilium_constructor(X0, L, NV, K, B, kap, bap, orient, ctype, nplus, nminus) result(self)
       real(real64), intent(in) :: X0(2) ! Length of cilium
       real(real64), intent(in) :: L ! Length of cilium
       integer(int32), intent(in) :: NV       
       real(real64), intent(in) :: K
       real(real64), intent(in) :: B
       real(real64), intent(in) :: kap
       real(real64), intent(in) :: bap
       logical, intent(in) :: orient
       logical, intent(in) :: ctype
       real(real64), intent(in) :: nplus
       real(real64), intent(in) :: nminus

       integer(int32) :: i,j,xid
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
       self%orient = orient
       self%ctype = ctype
       
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
                     self%Umid(self%NE,2), self%mdenmid(self%NE), &
                     self%SD(self%NE), self%SV(self%NE), self%motorf(self%NE), &
                     self%nplus(self%NE), self%nminus(self%NE))

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
       self%release = .false.

       ! Cilia motility
       self%SD = 0.0d0
       self%SV = 0.0d0
       self%motorf = 0.0d0

       ! Initialize bound motor fractions
       self%nplus = nplus
       self%nminus = nminus

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

    ! if (orient) then
    !    self%XV(:,1) = X0(1)
    !    self%XV(:,2) = X0(2) - self%dl * [(i, i=0, self%NE)]
    ! !    forall(xid=1:self%NV)
    !         ! self%XV(xid,1) = self%XV(self%NV-xid+1,1)
    !         ! self%XV(xid,2) = self%XV(self%NV-xid+1,2)
    ! !    end forall 
    ! end if

       ! Calculate the mid-points of the elements 
       self%XE = self%XV(1:self%NE,:) + (self%XV(2:self%NV,:) - self%XV(1:self%NE,:)) / 2.0d0
       self%XEOLD = self%XE
        ! Original anchor point location
        if (self%orient) then
            self%XM0 = self%XE(self%NE,:)
        else
            self%XM0 = self%XE(1,:)
        end if
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
        
        integer(int32) :: k, j
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

        ! Add moments to the tips
        if (present(mtip)) then
            ! mtip = self%B / (self%L * 2.0d0)
            ! mtip = self%B / self%rada
            ! self%mom(1)   = self%mom(1)  + mtip
            ! self%mom(self%NV)  = self%mom(self%NV) + mtip
            ! self%mden(self%NE) = self%mden(self%NE) + mtip
            self%mden(2) = self%mden(2) + mtip
            ! do j=2,self%NE
                ! self%mden(j) = self%mden(self%NE) + (mtip/self%NE)*(j*self%dl)
            ! end do
        end if
        ! Add the motor's contribution
        self%mden = self%mden + self%motorf
        ! print *, self%motorf

        ! Apply anchor point moment
        ! theta = PI/2.0d0 - atan2(self%that(1,2), self%that(1,1))
        ! !print *, "THe = ", theta
        
        ! momap = -self%bap * theta
        if (self%orient) then 
            ! theta = PI/2.0d0 - atan2(self%that(self%NE,2), self%that(self%NE,1))
            ! momap = -self%bap * theta
            ! self%mden(self%NE) = self%mden(self%NE) + (-momap)
            theta = PI/2.0d0 - atan2(self%that(self%NE,2), self%that(self%NE,1))
            momap = -self%bap * theta
            self%mden(self%NE) = self%mden(self%NE) + (-momap)
        else
            theta = PI/2.0d0 - atan2(self%that(1,2), self%that(1,1))
            momap = -self%bap * theta
            self%mden(1) = self%mden(1) + (-momap)
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
       
        ! fap(1) = self%kap * (self%XM0(1)-self%XE(1,1))
        ! fap(2) = self%kap * (self%XM0(2)-self%XE(1,2))
        !print *, "Fap = ", fap

        if (self%ctype) then ! Cilia is flow altering type (ctype = 1)
            if (self%orient) then ! Cilia is up side down (orient = 1) and flow altering type
                ! Depending upon the orientation the captive forces will change 
                fap(1) = self%kap * (self%XM0(1)-self%XE(self%NE,1)) ! Here we calculate the last node anchor forces
                fap(2) = self%kap * (self%XM0(2)-self%XE(self%NE,2)) ! Since the cilia is upside down    
                if (self%release) then ! Cilia is released
                    self%fden(self%NE,:)  = self%fden(self%NE,:) + fap ! Anchor point
                else ! Cilia is captive by horizontal springs along its nodes
                    self%fden(self%NE,:)  = self%fden(self%NE,:) + fap ! Anchor point
                    self%fden(1:self%NE-1,1)  = self%fden(1:self%NE-1,1) + self%kap * (self%XM0(1)-self%XE(1:self%NE-1,1))
                end if
            else ! Cilia is right side up and cilia is flow altering type
                fap(1) = self%kap * (self%XM0(1)-self%XE(1,1))
                fap(2) = self%kap * (self%XM0(2)-self%XE(1,2))
                if (self%release) then ! Cilia is released
                    self%fden(1,:)  = self%fden(1,:) + fap
                else ! Cilia is captive by horizontal springs along its nodes
                    self%fden(1,:)  = self%fden(1,:) + fap
                    self%fden(2:self%NE,1)  = self%fden(2:self%NE,1) + self%kap * (self%XM0(1)-self%XE(2:self%NE,1))
                end if
            end if

        else ! Cilia is flow detection type
            fap(1) = self%kap * (self%XM0(1)-self%XE(1,1)) ! Only apply anchor forces at the bottom
            fap(2) = self%kap * (self%XM0(2)-self%XE(1,2))
            self%fden(1,:)  = self%fden(1,:) + fap
        end if


        if (self%orient) then 
            fap(1) = self%kap * (self%XM0(1)-self%XE(self%NE,1))
            fap(2) = self%kap * (self%XM0(2)-self%XE(self%NE,2))

            if (self%release) then
                self%fden(self%NE,:)  = self%fden(self%NE,:) + fap ! Anchor point
            else
                self%fden(self%NE,:)  = self%fden(self%NE,:) + fap ! Anchor point
                self%fden(1:self%NE-1,1)  = self%fden(1:self%NE-1,1) + self%kap * (self%XM0(1)-self%XE(1:self%NE-1,1))
            end if
            
        else
            fap(1) = self%kap * (self%XM0(1)-self%XE(1,1))
            fap(2) = self%kap * (self%XM0(2)-self%XE(1,2))
            self%fden(1,:)  = self%fden(1,:) + fap
        end if

    end subroutine forces

    subroutine velocity(self)
        class(cilium), intent(inout) :: self

        self%U = ((self%F(2:self%NV,:) - self%F(1:self%NV-1,:))) / self%dl
    end subroutine velocity
    
    elemental subroutine update(self, dt)
        class(cilium), intent(inout) :: self
        real(real64), intent(in) :: dt
        
        self%XE(1:self%NE,1) = self%XE(1:self%NE,1) + dt * self%U(1:self%NE,1)
        self%XE(1:self%NE,2) = self%XE(1:self%NE,2) + dt * self%U(1:self%NE,2)

        ! self%phi(1:self%NE) = self%phi(1:self%NE) + dt * self%mden(1:self%NE)

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

    elemental impure subroutine calculate_motorf(self, dt, ma, mk, f0, fc, v0, mrho, mt, pi0, eps0)
        class(cilium), intent(inout) :: self
        real(real64), intent(in) :: dt, ma, mk, f0, fc, v0, mrho, mt, pi0, eps0

        real(real64) :: nbar, ntilde, phi0, pi01, eps01
        integer(int32) :: ii

        ! Clamped boundary condition
        phi0 = PI/2.0d0

        ! Calculate Sliding displacement, SD and Sliding velocity, SV
        self%SD = ma * (self%phi - phi0) 
        ! self%SD = ma * (self%phi - self%phi(1)) 
        self%SV = ma * self%mden 
        
        do ii = 1,self%NE
            nbar    = self%nplus(ii) + self%nminus(ii)
            ntilde  = self%nplus(ii) - self%nminus(ii)

            ! Calculate motor force
            self%motorf(ii) = f0*mrho*(nbar - (self%SV(ii)/v0) * ntilde) - mk*self%SD(ii)
        end do

        ! Multiply it with "-a" so that it can be added to the moment density term
        self%motorf = -ma * self%motorf
        ! self%motorf = ma * self%motorf
        ! self%motorf(self%NE) = 0.0d0 
        ! self%motorf(1:5) = 0.0d0 

        pi01  = pi0 
        eps01 = eps0
        ! Evolve the bound motor fractions
        do ii = 1,self%NE
            pi01  = 0.17 
            eps01 = 0.73
            self%nplus(ii) = self%nplus(ii) + dt*(pi01*(1 - self%nplus(ii)) - & 
                                            eps01*self%nplus(ii)*exp((f0/fc) *(1 - self%SV(ii)/v0)))
            
            pi01  = 0.25 
            eps01 = 0.75
            self%nminus(ii) = self%nminus(ii) + dt*(pi01*(1 - self%nminus(ii)) - & 
                                            eps01*self%nminus(ii)*exp((f0/fc) *(1 + self%SV(ii)/v0)))
        end do

        ! print *, self%nplus - self%nminus
        ! print *, self%motorf
    
    end subroutine calculate_motorf

end module mod_cilia    
