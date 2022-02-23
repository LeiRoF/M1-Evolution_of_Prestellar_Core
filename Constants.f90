    MODULE CONSTANTS

    IMPLICIT NONE

    INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(P=8)

!   Math constants

    REAL (KIND=dp), PUBLIC, PARAMETER :: pi  = 3.1415926535897932384_dp          ! Pi (!)
    REAL (KIND=dp), PUBLIC, PARAMETER :: sqpi = 1.77245385090551602_dp           ! Sqrt(pi)

!   Physical constants (CGS)
    
    REAL (KIND=dp), PUBLIC, PARAMETER :: c_lum  = 2.99791458e+10_dp              ! Velocity of light   [cm.s-1]
    REAL (KIND=dp), PUBLIC, PARAMETER :: h_planck = 6.62606876e-27_dp            ! Planck constant     [erg.s]
    REAL (KIND=dp), PUBLIC, PARAMETER :: kb = 1.3806503e-16_dp                   ! Boltzmann constant  [erg.K-1]
    REAL (KIND=dp), PUBLIC, PARAMETER :: NAvogadro = 6.022d23                    ! Avogadro number     [mol-1]
    REAL (KIND=dp), PUBLIC, PARAMETER :: mH = 1.67262178d-27                     ! Proton mass         [g]

    REAL (KIND=dp), PUBLIC, PARAMETER :: G = 6.67408d-8                          !  Gravitational constant [cm3.g-1.s-2]
    REAL (KIND=dp), PUBLIC, PARAMETER :: Msun = 2d33                             !  Mass of the sun [g]
    
    
!   Conversion factors
    REAL (KIND=dp), PUBLIC, PARAMETER :: G0 = 5.6e-14_dp                         ! G0
    REAL (KIND=dp), PUBLIC, PARAMETER :: eV_to_erg  = 1.602176487e-12_dp         ! eV->erg
    REAL (KIND=dp), PUBLIC, PARAMETER :: cm_to_A  = 1.0e+8_dp                    ! cm->Angström
    REAL (KIND=dp), PUBLIC, PARAMETER :: cmm1_to_eV  = 1.2398389635e-4_dp        ! cm-1->eV
    REAL (KIND=dp), PUBLIC, PARAMETER :: A_to_eV = h_planck &                    ! Angström->eV : careful to invert the variable ! (wl(eV)=A_to_eV/wl(A))
                                                * c_lum*cm_to_A &
                                                / eV_to_erg                      ! = 12398.377023296383
    REAL (KIND=dp), PUBLIC, PARAMETER :: parsec_to_cm = 3.08567758e18_dp         ! [cm.pc-1]

    REAL (KIND=dp), PUBLIC, PARAMETER :: DEGREE_TO_RADIAN = pi/180d0

!   Astro parameters
    REAL (KIND=dp), PUBLIC, PARAMETER :: mu = 1.4_dp                            ! Mean particle mass in unit of mH (proton mass)
    REAL (KIND=dp), PUBLIC, PARAMETER :: fractal = 0.3_dp                        ! Effective fraction of optical depth in Vis-UV RT, due to fractal structure of the ISM (see Varosi and Dwek 1999) 

    END MODULE CONSTANTS

