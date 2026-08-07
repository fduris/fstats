! =======================================================================
! MTRobust -- Native Theil-Sen and Siegel repeated-median robust linear
! fits, replacing
! R's `mblm` package (Median-Based Linear Models) used throughout
! `source_scripts/fitSpeed.r` (mblm(y~x), mblm(..., repeated=TRUE), and its
! `summary` read by `mblm2list`).
!
! ------------------------------------------------------------------------
! Pinned mblm semantics (read from the installed package, not assumed):
!   Rscript -e 'library(mblm); print(mblm); print(getS3method("summary","mblm"))'
!
! (a) DEFAULT of `repeated`: mblm's formal is `mblm(formula, dataframe,
!     repeated = TRUE)`, so the bare `mblm(y~x)` in fitSpeed.r uses
!     repeated=TRUE -> the SIEGEL repeated-median estimator. Other call
!     sites also pass repeated=TRUE explicitly. The app therefore only ever uses
!     Siegel; Theil-Sen (repeated=FALSE) is provided for completeness.
!
! (b) Coefficient definitions (mblm sorts by x first; because both estimators
!     reduce to order-independent medians of symmetric pairwise quantities,
!     the sort does not change the result, so it is omitted here):
!     * Theil-Sen (repeated=FALSE):
!         slope     = median over pairs i<j with x_j /= x_i of
!                     (y_j - y_i)/(x_j - x_i)
!         intercept = median over k of (y_k - slope*x_k)
!     * Siegel (repeated=TRUE):
!         for each point i: smedian_i = median over j /= i (x_j /= x_i) of
!                     (y_j - y_i)/(x_j - x_i)
!                           imedian_i = median over j /= i (x_j /= x_i) of
!                     (x_j*y_i - x_i*y_j)/(x_j - x_i)   [line through i,j]
!         slope     = median_i(smedian_i)
!         intercept = median_i(imedian_i)
!
! (c) summary.mblm coefficient SEs and sigma (what mblm2list reads):
!     * sigma = summary$sigma = sqrt( sum(residuals^2) / (n - 2) ), with
!       residuals = y - slope*x - intercept and rdf = df.residual = n - 2.
!       (This is the residual-standard-error form, NOT a MAD.)
!     * The coefficient-error column mblm2list reads (summary$coefficients[,2])
!       is the "MAD" column: c(mad(z$intercepts), mad(z$slopes)) where R's
!       mad(v) = 1.4826 * median(|v - median(v)|). The vectors z$slopes /
!       z$intercepts are the estimator's own working medians:
!         repeated=FALSE: z$slopes = all pairwise slopes,
!                         z$intercepts = y - slope*x (per point)
!         repeated=TRUE : z$slopes = smedian_i, z$intercepts = imedian_i
!       so seSlope = mad(<slope vector>), seIntercept = mad(<intercept
!       vector>). These are Wilcoxon/MAD statistics, not OLS standard errors.
!     (summary.mblm's V value / Pr(>|V|) columns are unused by mblm2list and
!     are not reproduced.)
!
! R's median: for m values sorted ascending, the middle one when m is odd and
! the mean of the two central ones when m is even -- replicated in MedianOf.
! ------------------------------------------------------------------------
! =======================================================================
module MTRobust

use, intrinsic :: iso_fortran_env, only : real64, int32, int64
use MTSelect, only : MedianInplace, MadInplace, MedianOf, MadOf
implicit none
private

public :: TheilSenFit, SiegelFit

! Largest sample the Theil-Sen estimator accepts, enforced by the
! ComputeTheilSenFit export. The estimator is defined over ALL n(n-1)/2 pairwise
! slopes and its MAD needs them materialised, so its memory is inherently
! quadratic: 4096 points is 8.4e6 pairs, a 67 MB buffer. That is the honest
! bound -- the general MAX_TRAIL_LENGTH of 100000 would be 5e9 pairs (40 GB) and
! could never run. Siegel (O(n) storage) is unaffected and keeps the general
! bound.
integer(int32), parameter, public :: MAX_THEILSEN_LENGTH = 4096

contains

! -----------------------------------------------------------------------
! ResidualSigma -- summary.mblm$sigma = sqrt( sum(residuals^2) / (n - 2) ).
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
! TheilSenFit -- mblm(y~x, repeated=FALSE). Preconditions: 3 <= n <=
! MAX_THEILSEN_LENGTH and
! not all x equal (guaranteed by the caller); the pairwise loop then yields
! >= 1 slope. `status` is 0 on success and 1 if the pairwise-slope buffer
! could not be allocated.
! -----------------------------------------------------------------------
subroutine TheilSenFit(n, x, y, slope, intercept, sigma, seSlope, seIntercept, &
        status)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n), y(n)
    real(real64), intent(out) :: slope, intercept, sigma, seSlope, seIntercept
    integer(int32), intent(out) :: status

    ! The pairwise-slope buffer is ALLOCATABLE, not automatic: at the documented
    ! maximum it is 67 MB, which would overflow any ordinary thread stack. The
    ! pair count is formed in int64 because n*(n-1) overflows a 32-bit integer
    ! above n = 65536 -- as an automatic array bound that silently produced a
    ! wrapped, meaningless extent (n = 70000 asked for 302,481,352 elements
    ! instead of 2,449,965,000) and, for wraps that land negative, a zero-size
    ! array written far out of bounds.
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

    ! MAD columns read by mblm2list: on the pairwise slopes and the per-point
    ! intercepts respectively (see the semantics block above). mad's own centre is
    ! median(<vector>), i.e. exactly the estimate just computed, so `slopes` can be
    ! consumed in place rather than copied a second time.
    call MadInplace(slopes(1:ns), ns, slope, seSlope)
    seIntercept = MadOf(intercepts, n)

    deallocate(slopes, intercepts)
end subroutine TheilSenFit

! -----------------------------------------------------------------------
! SiegelFit -- mblm(y~x, repeated=TRUE). Precondition: n >= 3 and not all
! x equal
! (guaranteed by the caller); every point then has >= 1 partner with
! distinct x, so each per-point median is defined. Storage is O(n) -- only
! the per-point medians are materialised, never the full pair set -- so this
! estimator keeps the general MAX_TRAIL_LENGTH bound. `status` is 0 on
! success and 1 if the working buffers could not be allocated.
! -----------------------------------------------------------------------
subroutine SiegelFit(n, x, y, slope, intercept, sigma, seSlope, seIntercept, &
        status)
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

    ! MAD columns read by mblm2list: on the per-point median vectors. Their mad
    ! centres are the estimates just computed, so both vectors are consumed in
    ! place.
    call MadInplace(smedians, n, slope, seSlope)
    call MadInplace(imedians, n, intercept, seIntercept)

    deallocate(smedians, imedians, pslopes, pints)
end subroutine SiegelFit

end module MTRobust
