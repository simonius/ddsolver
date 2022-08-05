module consts
        implicit none
! Fortran is not case sensitive
        Character(LEN=20) :: filename = "hdf5test.h5"
        Character(LEN=20) :: udsetname = "u"
        Character(LEN=20) :: alpdsetname = "alpha"
        Character(LEN=20) :: edsetname = "dE"

        integer, parameter :: sols = 20
        integer, parameter :: ntest = 200
        integer, parameter :: fac = 1
        integer, parameter :: fdr = 10
        integer, parameter :: facp = 10
        integer, parameter :: n = fac*ntest
        integer, parameter :: ncons = 3
        integer, parameter :: ktest = fdr*ntest
        integer, parameter :: k = fdr*n
        double precision, parameter :: xmin = 0.0
        double precision, parameter :: xmax = 10.0
        double precision, parameter :: CFL = 0.025
        double precision, parameter :: lambda = CFL
        double precision, parameter :: dx = (xmax-xmin) / n
        double precision, parameter :: dt = lambda * dx
        double precision, parameter :: gamma = 1.4
        double precision, parameter :: PI=4.D0*DATAN(1.D0)
        double precision, parameter :: tend = k*dt
        logical, parameter :: diag = .False.


        double precision, parameter :: rhomax = 5.0D0
        double precision, parameter :: vmax = 5.0D0
        double precision, parameter :: pmax = 5.0D0
        double precision, parameter :: GTeps = 1.0D-14
        double precision, parameter :: rhoeps = 0.05D0
        double precision, parameter :: peps = 0.05D0
        double precision, parameter :: ppeps = 1.0D-5

        contains

        function warp(i, imax)
                integer :: i, imax
                integer :: warp
                warp = MODULO(i-1, imax)+1
        end function warp
end module consts
