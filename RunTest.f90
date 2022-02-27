PROGRAM TEST
   USE Constants
   USE PRESTEL
   USE UnitaryTests
   IMPLICIT NONE

   INTEGER, parameter             :: N =100
   REAL (KIND=dp)                 :: Rmin = 1.0_dp, Rmax = 1.0e3_dp
   real(kind=dp)                  :: Temp = 10 ! K

   integer                        :: i
   real (kind=dp)                 :: cs, dr, dt, tff, mass, res, sigma, rho_prestel, volume
   real (kind=dp), dimension(N)   :: r, distrib, t, buffer
   real (kind=dp), dimension(N,N) :: rho, v

   ! Test of all subroutines
   
   call EXECUTE_COMMAND_LINE('mkdir results')
   call EXECUTE_COMMAND_LINE('mkdir "results/unitary_test"')
   call RunAllUnitaryTests()

   ! Analysis

   mass = Msun
   volume = 4/3 * pi * (Rmax**3 - Rmin**3)
   rho_prestel = mass/volume
   call time_free_fall(rho_prestel, tff)
   call get_cs(Temp, cs)
   call linspace(Rmin, Rmax, N, r, dr)
   call linspace(0e0_dp,tff, N, t, dt)

   sigma = N / 2
   do i = 1, N
      distrib(i) = 1/(sqrt(2*pi) * sigma) * exp(-(i-1)**2/(2*sigma**2)) ! Gaussian density distribution
      v(1,:)   = 0                                                            ! No initial speed
   end do

   !call density_1D(distrib, mass, N, r, dr, buffer)                           ! Normalize the distribution to conserve the total mass
   rho(1,:) = distrib !buffer

   open (unit = 40, file = "results/space.dat")
   write(40,*) r
   close(40)

   open (unit = 40, file = "results/time.dat")
   write(40,*) t
   close(40)

   open (unit = 40, file = "results/rho.dat")
   write(40,*) rho(1,:)
   open (unit = 41, file = "results/v.dat")
   write(41,*) v(1,:)

   do i=2,N
      call next_speed(v(i-1,:), rho(i-1,:), r, N, dr, dt, cs, v(i,:))
      call next_density(rho(i-1,:), v(i-1,:), r, N, dr, dt, rho(i,:))
      write(40,*) rho(i,:)
      write(41,*) v(i,:)
   end do

   close(40)
   close(41)

   print *, " "
   print *, "-------------------------------"
   print *, "   PRESTELLAR GAZ COLLAPSING   "
   print *, "-------------------------------"
   print *, " "
   print *, "Mass             = ", mass / Msun, " time the mass of the sun, "
   print *, "Radius           = ", Rmax / 100000 ," km (=", Rmax / parsec_to_cm, " parsec), "
   print *, "Free fall time   = ", tff / 3600 / 24 / 365.25, " years."
   print *, "Density          = ", rho_prestel, " g.cm-3"
   print *, "Iso. sound speed = ", cs


END PROGRAM TEST
