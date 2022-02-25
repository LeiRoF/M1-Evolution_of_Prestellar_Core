PROGRAM TEST

   USE Constants
   USE PRESTEL

   IMPLICIT NONE

   INTEGER, parameter :: N=100
   REAL (KIND=dp) :: Rmin=1.0_dp, Rmax=1.0e3_dp

   integer :: i=0, tmp
   real (kind=dp) :: cs, Temp = 10, dr, dt = 1, tff, rho_m3D = 1, sigma, mass, buffer2, res, avd1D
   real (kind=dp), dimension(N) :: r, distrib, t, buffer
   real (kind=dp), dimension(N,N) :: rho, v
   character(len=1024) :: str_i, filename
   
   print *, " "
   print *, "-------------------------------"
   print *, "   PRESTELLAR GAZ COLLAPSING   "
   print *, "-------------------------------"

   !Rmax = 0.1/2 * parsec_to_cm ! [cm]
   mass = Msun ! [g]
   rho_m3D = mass / (4./3 * pi * Rmax**3)  ! [g.cm-3]


   call time_free_fall(rho_m3D, tff)
   call get_cs(Temp, cs)
   
   call linspace(Rmin, Rmax, N, r, dr)
   call linspace(0e0_dp,tff, N, t, dt)

   print *, " "
   print *, "Mass             = ", mass / Msun, " time the mass of the sun, "
   print *, "Radius           = ", Rmax / 100000 ," km (=", Rmax / parsec_to_cm, " parsec), "
   print *, "Free fall time   = ", tff / 3600 / 24 / 365.25, " years."
   print *, "Density          = ", rho_m3D, " g.cm-3"
   print *, "Iso. sound speed = ", cs



   ! __________________________________________________
   ! Distribution of mass in the gas cloud

   sigma = N / 2
   do i = 1, N
      !distrib(i) = 1/(sqrt(2*pi) * sigma) * exp(-((i-1)-N/2)**2/(2*sigma**2)) ! Gaussian distribution
      distrib(i) = 1
   end do

   !call integrate_on_sphere(D,N,r,dr,avd1D)

   call density_1D(distrib, mass, N, r, dr, rho(1,:))

   ! __________________________________________________
   ! Initial vectors
   
   do i = 1,N
      v(1,i) = 0
      rho(1,i) = 1 !rho_m3D
   end do


   ! __________________________________________________
   ! Computation

   call total_mass(rho(1,:), N, r, dr, mass)

   print *, " "
   print *, "mass = ", mass/Msun, "*Msun // TO DO (la masse devrait etre proche de 1 masse solaire)"

   call EXECUTE_COMMAND_LINE("mkdir results")
   
   call derive(rho(1,:), dr, N, buffer)
   open (unit = 40, file = "results/d_rho.dat")
   write(40,*) buffer
   close(40)

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

   !v(1,:) = r
   do i=2,N
      call next_speed(v(i-1,:), rho(i-1,:), r, N, dr, dt, cs, v(i,:))
      call next_density(rho(i-1,:), v(i-1,:), r, N, dr, dt, rho(i,:))
      write(40,*) rho(i,:)
      write(41,*) v(i,:)
   end do

   close(40)
   close(41)
   
   print *, " "
END PROGRAM TEST
