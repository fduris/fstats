! =======================================================================
! MTRobust -- Theil-Sen and Siegel repeated-median robust linear fits.
!
! Both estimate a straight line y = intercept + slope*x from medians of
! pairwise quantities, so a minority of grossly wrong points cannot drag
! the fit. Only pairs with distinct abscissae contribute.
!
!   Theil-Sen
!     slope     = median over pairs i<j of (y_j - y_i)/(x_j - x_i)
!     intercept = median over k of (y_k - slope*x_k)
!
!   Siegel (repeated median, the estimator the speed fit uses)
!     per point i, over partners j /= i:
!       smedian_i = median of (y_j - y_i)/(x_j - x_i)
!       imedian_i = median of (x_j*y_i - x_i*y_j)/(x_j - x_i)
!     slope     = median_i(smedian_i)
!     intercept = median_i(imedian_i)
!
! Both report sigma = sqrt(sum(residuals^2)/(n-2)) and, as the per-
! coefficient error, the MAD of the estimator's own working vector: the
! pairwise slopes and per-point intercepts for Theil-Sen, the two
! per-point median vectors for Siegel. Those are MAD statistics, not OLS
! standard errors.
!
! Stateless: no SAVE, no module variables, all working storage local or
! caller-passed.
! =======================================================================
module MTRobust

use, intrinsic :: iso_fortran_env, only : real64, int32, int64
use MTSelect, only : MedianInplace, MadInplace, MedianOf, MadOf
implicit none
private

public :: TheilSenFit, SiegelFit

! Largest sample the Theil-Sen estimator accepts, enforced by the
! ComputeTheilSenFit export. All n(n-1)/2 pairwise slopes must be
! materialised for the MAD, so storage is quadratic: 4096 points is 8.4e6
! pairs, a 67 MB buffer. Siegel needs only O(n) and keeps the general
! trail bound.
integer(int32), parameter, public :: MAX_THEILSEN_LENGTH = 4096

contains

! -----------------------------------------------------------------------
! ResidualSigma -- sqrt( sum(residuals^2) / (n - 2) ).
! Precondition: n >= 3, so the denominator is positive.
! -----------------------------------------------------------------------
pure function ResidualSigma(n, x, y, slope, intercept) result(sigma)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n), y(n), slope, intercept
    real(real64) :: sigma

    real(real64) :: r(n)

    r = y - slope * x - intercept
    sigma = sqrt(sum(r * r) / real(n - 2, real64))
end function ResidualSigma

! -----------------------------------------------------------------------
! TheilSenFit -- Theil-Sen fit, defined in the module header.
!
! Preconditions, guaranteed by the caller: 3 <= n <= MAX_THEILSEN_LENGTH
! and not all x equal, so the pairwise loop yields at least one slope.
! `status` is 0 on success, 1 if the pairwise-slope buffer could not be
! allocated.
! -----------------------------------------------------------------------
subroutine TheilSenFit(n, x, y, slope, intercept, sigma, seSlope, seIntercept, status)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n), y(n)
    real(real64), intent(out) :: slope, intercept, sigma, seSlope, seIntercept
    integer(int32), intent(out) :: status

    ! Allocatable, not automatic: at the maximum n this is 67 MB, which would
    ! overflow an ordinary thread stack. The pair count is int64 because
    ! n*(n-1) overflows int32 above n = 65536.
    real(real64), allocatable :: slopes(:)
    real(real64), allocatable :: intercepts(:)
    integer(int64) :: npair
    integer(int32) :: i, j, ns, ierr

    slope = 0.0_real64
    intercept = 0.0_real64
    sigma = 0.0_real64
    seSlope = 0.0_real64
    seIntercept = 0.0_real64
    status = 0

    npair = int(n, int64) * int(n - 1, int64) / 2_int64
    allocate(slopes(npair), intercepts(n), stat=ierr)
    if (ierr /= 0) then
        status = 1
        return
    end if

    ! All pairwise slopes over i < j with x_j /= x_i.
    ns = 0
    do i = 1, n - 1
        do j = i + 1, n
            if (x(j) /= x(i)) then
                ns = ns + 1
                slopes(ns) = (y(j) - y(i)) / (x(j) - x(i))
            end if
        end do
    end do
    call MedianInplace(slopes(1:ns), ns, slope)

    ! Per-point intercepts about the single fitted slope.
    intercepts = y - slope * x
    intercept = MedianOf(intercepts, n)

    sigma = ResidualSigma(n, x, y, slope, intercept)

    ! Each MAD centres on the estimate just computed, so both working vectors
    ! are consumed in place rather than copied.
    call MadInplace(slopes(1:ns), ns, slope, seSlope)
    seIntercept = MadOf(intercepts, n)

    deallocate(slopes, intercepts)
end subroutine TheilSenFit

! -----------------------------------------------------------------------
! SiegelFit -- Siegel repeated-median fit, defined in the module header.
!
! Preconditions, guaranteed by the caller: n >= 3 and not all x equal, so
! every point has at least one partner with distinct x and each per-point
! median is defined. Only the per-point medians are materialised, never
! the full pair set, so storage is O(n). `status` is 0 on success, 1 if
! the working buffers could not be allocated.
! -----------------------------------------------------------------------
subroutine SiegelFit(n, x, y, slope, intercept, sigma, seSlope, seIntercept, status)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n), y(n)
    real(real64), intent(out) :: slope, intercept, sigma, seSlope, seIntercept
    integer(int32), intent(out) :: status

    ! Allocatable rather than automatic so nothing is sized onto the caller's
    ! stack from an n it does not control.
    real(real64), allocatable :: smedians(:), imedians(:)
    real(real64), allocatable :: pslopes(:), pints(:)   ! per-point (<= n-1 used)
    integer(int32) :: i, j, k, ierr

    slope = 0.0_real64
    intercept = 0.0_real64
    sigma = 0.0_real64
    seSlope = 0.0_real64
    seIntercept = 0.0_real64
    status = 0

    allocate(smedians(n), imedians(n), pslopes(n), pints(n), stat=ierr)
    if (ierr /= 0) then
        status = 1
        return
    end if

    do i = 1, n
        k = 0
        do j = 1, n
            if (x(j) /= x(i)) then
                k = k + 1
                pslopes(k) = (y(j) - y(i)) / (x(j) - x(i))
                pints(k) = (x(j) * y(i) - x(i) * y(j)) / (x(j) - x(i))
            end if
        end do
        smedians(i) = MedianOf(pslopes(1:k), k)
        imedians(i) = MedianOf(pints(1:k), k)
    end do
    slope = MedianOf(smedians, n)
    intercept = MedianOf(imedians, n)

    sigma = ResidualSigma(n, x, y, slope, intercept)

    ! Each MAD centres on the estimate just computed, so both vectors are
    ! consumed in place.
    call MadInplace(smedians, n, slope, seSlope)
    call MadInplace(imedians, n, intercept, seIntercept)

    deallocate(smedians, imedians, pslopes, pints)
end subroutine SiegelFit

end module MTRobust
