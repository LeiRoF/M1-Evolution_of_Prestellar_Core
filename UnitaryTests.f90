module UnitaryTests
    use Constants
    use Prestel

    implicit none

    CONTAINS

    subroutine RunAllUnitaryTests()
        
        call test_density_1D()
        call test_density_to_nH()
        call test_derive()
        call test_get_cs()
        call test_grav_force()
        call test_integrate_on_sphere()
        call test_linspace()
        call test_mass_accretion_rate()
        call test_next_density()
        call test_next_speed()

    end subroutine RunAllUnitaryTests

    subroutine test_density_1D()

        integer, parameter          :: N = 100
        real(kind=dp)               :: mass = 1.0, dr, sigma
        real(kind=dp), dimension(N) :: distrib, r, rho

        call linspace(0.0_dp, 1.0_dp, N, r, dr)

        sigma = 1.0_dp / 2
        distrib = 1/(sqrt(2*pi) * sigma) * exp(-r**2/(2*sigma**2)) ! Gaussian distribution

        call density_1D(distrib, mass, N, r, dr, rho)

        open (unit = 40, file = "results/unitary_test/density_1D.dat")
        write(40,*) rho
        close(40)

    end subroutine test_density_1D

    subroutine test_density_to_nH()
        real(kind=dp) :: rho = 3.0, res

        call density_to_nH(rho,res)

        open (unit = 40, file = "results/unitary_test/density_to_nH.dat")
        write(40,*) res
        close(40)

    end subroutine test_density_to_nH

    subroutine test_derive()
        integer, parameter          :: N = 100
        real(kind=dp)               :: dx
        real(kind=dp), dimension(N) :: x, res, F

        call linspace(0.0_dp, 1.0_dp * 5 * pi, N, x, dx)

        F = sin(x) / x

        call derive(F, dx, N, res)

        open (unit = 40, file = "results/unitary_test/derive.dat")
        write(40,*) res
        close(40)

    end subroutine test_derive

    subroutine test_get_cs()
        real(kind=dp) :: T = 10.0, res

        call get_cs(T, res)

        open (unit = 40, file = "results/unitary_test/get_cs.dat")
        write(40,*) res
        close(40)

    end subroutine test_get_cs

    subroutine test_grav_force()
        integer, parameter          :: N = 100
        real(kind=dp)               :: dr, dphi
        real(kind=dp), dimension(N) :: rho, r

        call linspace(0.0_dp, 1.0_dp, N, r, dr)
        rho = exp(-r)

        call grav_force(rho, N, r, dr, dphi)

        open (unit = 40, file = "results/unitary_test/grav_force.dat")
        write(40,*) dphi
        close(40)

    end subroutine test_grav_force

    subroutine test_integrate_on_sphere()
        integer, parameter          :: N = 100
        real(kind=dp)               :: dr, res
        real(kind=dp), dimension(N) :: D, r

        call linspace(0.0_dp, 5.0_dp, N, r, dr)

        D = exp(-r)

        call integrate_on_sphere(D, N, r, dr, res)

        open (unit = 40, file = "results/unitary_test/integrate_on_sphere.dat")
        write(40,*) res
        close(40)

    end subroutine test_integrate_on_sphere

    subroutine test_linspace()
        integer, parameter          :: N = 100
        real(kind=dp)               :: dr
        real(kind=dp), dimension(N) :: r

        call linspace(0.0_dp, 5.0_dp, N, r, dr)

        open (unit = 40, file = "results/unitary_test/linspace.dat")
        write(40,*) r
        close(40)

        open (unit = 40, file = "results/unitary_test/linspace_step.dat")
        write(40,*) dr
        close(40)

    end subroutine test_linspace

    subroutine test_mass_accretion_rate()
        integer, parameter            :: N = 100
        integer                       :: i
        real(kind=dp)                 :: dr, dt, y
        real(kind=dp), dimension(N)   :: r, t , res, x
        real(kind=dp), dimension(N,N) :: rho

        call linspace(0.0_dp, 1.0_dp, N, r, dr)
        call linspace(0.0_dp, 1.0_dp, N, t, dt)

        do i=1,N
            x = (r-0.5)*10
            y = (t(i)-0.5)*10
            rho(i,:) = (x**2 + y - 11)**2 + (x + y**2 -7)**2
        end do

        call mass_accretion_rate(rho, 20, dt, N, N, res)

        open (unit = 40, file = "results/unitary_test/mass_accretion_rate_rho.dat")
        write(40,*) rho
        close(40)

        open (unit = 40, file = "results/unitary_test/mass_accretion_rate.dat")
        write(40,*) res
        close(40)

    end subroutine test_mass_accretion_rate


    subroutine test_next_density()
        integer, parameter            :: N = 100
        integer                       :: i
        real(kind=dp)                 :: dr, dt=1.0
        real(kind=dp), dimension(N)   :: r, res, rho, v

        call linspace(0.0_dp, 1.0_dp, N, r, dr)
        rho = exp(-r)
        v = sin(r*4*pi)

        call next_density(rho, v, r, N, dr, dt, res)

        open (unit = 40, file = "results/unitary_test/next_density.dat")
        write(40,*) res
        close(40)

        open (unit = 40, file = "results/unitary_test/next_density_v.dat")
        write(40,*) v
        close(40)

        open (unit = 40, file = "results/unitary_test/next_density_rho.dat")
        write(40,*) rho
        close(40)
    end subroutine test_next_density

    subroutine test_next_speed()
        integer, parameter            :: N = 100
        integer                       :: i
        real(kind=dp)                 :: dr, dt=1.0, cs = 3.0
        real(kind=dp), dimension(N)   :: r, res, rho, v

        call linspace(0.0_dp, 1.0_dp, N, r, dr)
        rho = exp(-r)
        v = sin(r*4*pi)

        call next_speed(v, rho, r, N, dr, dt, cs, res)
        
        open (unit = 40, file = "results/unitary_test/next_speed.dat")
        write(40,*) res
        close(40)
    end subroutine test_next_speed


end module UnitaryTests