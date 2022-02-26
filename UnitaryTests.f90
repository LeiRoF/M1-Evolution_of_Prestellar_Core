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

end module UnitaryTests