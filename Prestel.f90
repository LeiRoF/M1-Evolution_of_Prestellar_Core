!
!===========================
   MODULE Prestel
!===========================

   USE Constants

   IMPLICIT NONE

   CONTAINS

!===============================================================================
   subroutine mean_density_1D(shape, mass, N, r, dr, rho) ! OK.
!===============================================================================
      ! Get gravitational force

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: shape ! shape of the density distribution
      real (kind=dp)               :: mass  ! mass of the system
      integer                      :: N     ! array dimension
      real (kind=dp), dimension(N) :: r     ! space
      real (kind=dp)               :: dr    ! space step

      ! Outputs
      real (kind=dp), dimension(N) :: rho   ! distribution of the density 

      ! Private
      real (kind=dp)               :: tmp   ! buffer

      ! __________________________________________________
      ! Instructions

      call integrate_on_sphere(shape, N, r, dr, tmp)

      rho = shape / tmp * mass

   end subroutine mean_density_1D

!===============================================================================
   subroutine integrate_on_sphere(D, N, r, dr, res) ! OK.
!===============================================================================
      ! Integrate a 1D distribution on a sphere

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: D    ! 1D distribution array to integrate
      integer                      :: N    ! array dimension
      real (kind=dp), dimension(N) :: r    ! space array
      real (kind=dp)               :: dr   ! space step

      ! Outputs
      real (kind=dp)               :: res ! result

      ! Private
      integer                      :: i    ! index

      ! __________________________________________________
      ! Instructions

      res = 4 * pi * sum(D * r**2 * dr)  ! integral on a sphere

   end subroutine integrate_on_sphere

!===============================================================================
   subroutine total_mass(rho, N, r, dr, res) ! OK.
!===============================================================================
      ! Get gravitational force

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: rho  ! density array at a given time
      integer                      :: N    ! array dimension
      real (kind=dp), dimension(N) :: r    ! space array
      real (kind=dp)               :: dr   ! space step

      ! Outputs
      real (kind=dp)               :: res  ! total mass considered

      ! Private
      integer                      :: i    ! index

      ! __________________________________________________
      ! Instructions

      call integrate_on_sphere(rho,N,r,dr,res)

   end subroutine total_mass

!===============================================================================
   subroutine grav_force(rho, N, r, dr, dphi)
!===============================================================================
      ! Get gravitational force

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: rho  ! density array
      real (kind=dp), dimension(N) :: r    ! space array
      integer                      :: N    ! array dimension
      real (kind=dp)               :: dr   ! space step

      ! Outputs
      real (kind=dp)               :: dphi ! gravitational force

      ! Private
      real (kind=dp)               :: tmp  ! buffer

      ! __________________________________________________
      ! Instructions

      !dphi = 4 * pi * G * sum(rho) * dr ! 1D

      call integrate_on_sphere(rho, N, r, dr, tmp)
      dphi = 4 * pi * G * tmp ! 3D

   end subroutine grav_force

!===============================================================================
   subroutine get_cs(T, res) ! OK.
!===============================================================================
      ! Get isothermal speed of sound

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp) :: T   ! temprature [K]

      ! Outputs
      real (kind=dp) :: res ! isothermal speed of sound

      ! __________________________________________________
      ! Instructions

      res = sqrt(kb * T / (mu * mH))

   end subroutine get_cs

!===============================================================================
   subroutine linspace(rmin, rmax, N, res, dr) ! OK.
!===============================================================================
      ! Get a space representation in an array

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      integer                      :: N    ! number of discrete space points
      real (kind=dp)               :: rmin ! starting position
      real (kind=dp)               :: rmax ! stop position

      ! Outputs
      real (kind=dp), dimension(N) :: res  ! space array
      real (kind=dp)               :: dr   ! space step

      ! Private
      integer                      :: i    ! index

      ! __________________________________________________
      ! Instructions

      dr = (rmax - rmin)/(N-1)

      do i=1,N
         res(i) = rmin + (i-1)*dr
      end do

   end subroutine linspace

!===============================================================================
   subroutine next_state(rho, v, r, N, dr, dt, cs, new_rho, new_v)
!===============================================================================
      ! Simulation evolution of the state at time i (get the state at time i+1)

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: rho     ! array containing densities at time i
      real (kind=dp), dimension(N) :: v       ! array containing speeds    at time i
      real (kind=dp), dimension(N) :: r       ! array representing space
      integer                      :: N       ! dimension of the arrays
      real (kind=dp)               :: cs      ! isothermal speed of sound
      real (kind=dp)               :: dt      ! time step
      real (kind=dp)               :: dr      ! space step

      ! Outputs
      real (kind=dp), dimension(N) :: new_rho ! returning density          at time i+1
      real (kind=dp), dimension(N) :: new_v   ! returning speeds           at time i+1

      ! __________________________________________________
      ! Instructions

      call next_rho(rho, v, r, N, dr, dt, new_rho)
      call next_v(v, rho, r, N, dr, dt, cs, new_v)

   end subroutine next_state

!===============================================================================
   subroutine next_rho(rho, v, r, N, dr, dt, res)
!===============================================================================
      ! Get the density array at time i+1

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: rho ! array containing densities at time i
      real (kind=dp), dimension(N) :: v   ! array containing speeds    at time i
      real (kind=dp), dimension(N) :: r   ! array representing space
      integer                      :: N   ! dimension of the arrays
      real (kind=dp)               :: dt  ! time step
      real (kind=dp)               :: dr  ! space step

      ! Outputs
      real (kind=dp), dimension(N) :: res ! returning array rho        at time i+1

      ! Private
      real (kind=dp), dimension(N) :: tmp ! buffer

      ! __________________________________________________
      ! Instructions

      tmp = r**2 * rho * v
      call derive(tmp, dr, N, tmp)
      
      res = rho - ( tmp / r**2 ) * dt

   end subroutine next_rho

!===============================================================================
   subroutine next_v(v, rho, r, N, dr, dt, cs, res)
!===============================================================================
      ! Get the speed array at time i+1

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: v   ! array containing speeds    at time i
      real (kind=dp), dimension(N) :: rho ! array containing densities at time i
      real (kind=dp), dimension(N) :: r   ! array representing space
      integer                      :: N   ! dimension of the arrays
      real (kind=dp)               :: cs  ! isothermal speed of sound
      real (kind=dp)               :: dt  ! time step
      real (kind=dp)               :: dr  ! space step

      ! Outputs
      real (kind=dp), dimension(N) :: res ! returning speeds           at time i+1

      ! Private
      real (kind=dp), dimension(N) :: tmp1, tmp2 ! buffers
      real (kind=dp)               :: tmp3       ! buffers
      integer                      :: i          ! index

      ! __________________________________________________
      ! Instructions

      call derive(v, dr, N, tmp1)
      call derive(rho, dr, N, tmp2)
      
      call grav_force(rho,N,r,dr,tmp3)

      res = v + (-v * tmp1 - (cs**2/rho) * tmp2 - tmp3) * dt

   end subroutine next_v

!===============================================================================
   subroutine derive_at(f, x, dx, N, res) ! OK.
!===============================================================================
      ! Get the derivative at a given point using an arry containing some discret values of a function.

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: f     ! array
      integer                      :: x     ! position in the array where you want to get the derivative
      real (kind=dp)               :: dx    ! step of the array
      integer                      :: N     ! size of the array

      ! Outputs
      real (kind=dp)               :: res   ! value of the derivative

      ! __________________________________________________
      ! Instructions
      
      res = 0

      ! Not at the begining of the array 
      if (x > 1) then
         res = res + (f(x) - f(x-1)) / dx
      end if

      ! Not at the end of the array
      if (x < N) then
         res = res + (f(x+1) - f(x)) / dx
      end if

      ! Not at extremums (centred finite derivative method)
      if (x > 1 .and. x < N) then
         res = res / 2
      end if

   end subroutine derive_at

!===============================================================================
   subroutine derive(f, dx, N, res) ! OK.
!===============================================================================
      ! Get the derivated function. A function (and it's derivated form) is represented by an array containing discret values of this function.

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp), dimension(N) :: f     ! array
      real (kind=dp)               :: dx    ! step of the array
      integer                      :: N     ! size of the array

      ! Ouputs
      real (kind=dp), dimension(N) :: res   ! derivated array

      ! Private
      real (kind=dp)               :: tmp   ! buffer
      integer                      :: i     ! index

      ! __________________________________________________
      ! Instructions

      do i=1,N
         call derive_at(f,i,dx,N,tmp)
         res(i) = tmp
      end do

   end subroutine derive

!===============================================================================
   subroutine time_free_fall(rho_0, t) ! OK.
!===============================================================================
      ! Get the time of free fall for a gas that have an initial density of rho_0

      implicit none

      ! __________________________________________________
      ! Declarations

      ! Inputs
      real (kind=dp) :: rho_0 ! initial density [g.cm-3]
      
      ! Ouputs
      real (kind=dp) :: t     ! time of free fall [s]
      
      ! __________________________________________________
      ! Instructions

      t = sqrt(3 * pi / (32 * G * rho_0))

   end subroutine time_free_fall

! ____________________________________________________________________________________________________
   END MODULE Prestel
