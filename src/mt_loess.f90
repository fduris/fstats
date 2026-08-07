! =======================================================================
! MTLoess -- Locally-weighted polynomial regression (loess), evaluating an
! exact local fit at every input abscissa.
!
! Per evaluation abscissa x_i, with no robustness iterations:
!   q     = min(n, floor(span*n + Q_TOL))   neighbourhood size
!   dmax  = the q-th smallest |x_j - x_i|   bandwidth
!   w_j   = (1 - (|x_j-x_i|/dmax)^3)^3 for |x_j-x_i| < dmax, else 0
! then a weighted polynomial of the requested degree is fitted in the
! CENTRED abscissa (x_j - x_i) by weighted least squares. Centring makes
! the intercept the fitted value AT x_i, so no re-evaluation is needed.
! When q = n the farthest neighbour sits at d = dmax and gets weight
! exactly 0 via (1 - 1^3)^3.
!
! Fitting each point directly, rather than interpolating a surface over
! kd-tree cells, is what keeps this stateless: no SAVE, no module state,
! no I/O, all working storage local to the call.
! =======================================================================
module MTLoess

use, intrinsic :: iso_fortran_env, only : real64, int32
use MTSelect, only : SelectKth
implicit none
private

public :: LoessFit

! Rounding tolerance added before truncating span*n to the neighbourhood
! size. Not cosmetic: span*n is generally inexact in binary, so a bare
! floor drops a whole neighbour whenever the product lands a hair below an
! integer -- span 0.35 with n = 180 gives 62.99999999999999, hence q = 62
! where 63 is intended. That one neighbour moves the fitted values by
! 1.6e-3 * sd(y). Ordinary (span, n) pairs are affected, so the tolerance
! is required for correctness.
real(real64), parameter :: Q_TOL = 1.0e-5_real64

contains

! -----------------------------------------------------------------------
! LoessFit -- Direct-surface loess of `y` on `x` at the input abscissae.
! `degree` is 1 or 2; `span` in (0, 1] selects the neighbourhood
! fraction, a span above 1 being clamped to the full sample by the q =
! min(n, ...) below. On return `ys` holds the fitted values and `status`
! is 0 (computed) or 2 (a local neighbourhood was degenerate: fewer than
! degree+1 positive-weight points, all-coincident abscissae in the
! neighbourhood, or a singular weighted normal matrix). Such a
! neighbourhood surfaces status=2 rather than an unreliable smoothed
! value. Preconditions (n >= 1, degree in {1,2}, span > 0) are enforced
! by the ComputeLoessFit export.
! -----------------------------------------------------------------------
subroutine LoessFit(n, x, y, span, degree, ys, status)
    external :: dposv   ! LAPACK SPD solve (resolved via the linked liblapack)
    integer(int32), intent(in) :: n, degree
    real(real64), intent(in) :: x(n), y(n), span
    real(real64), intent(out) :: ys(n)
    integer(int32), intent(out) :: status

    real(real64), allocatable :: d(:), dsel(:)
    integer(int32) :: p, q, i, j, k, kk, npos, info
    real(real64) :: dmax, dx, wj
    real(real64) :: xtwx(degree + 1, degree + 1), xtwy(degree + 1), basis(degree + 1)

    status = 0
    ys = 0.0_real64
    p = degree + 1

    ! Neighbourhood size q = min(n, floor(span*n + Q_TOL)); q >= p is needed for
    ! a determined local polynomial (the per-point npos check below enforces it).
    q = min(n, int(floor(span * real(n, real64) + Q_TOL), int32))
    if (q < 1) then
        status = 2
        return
    end if

    allocate(d(n), dsel(n))

    do i = 1, n
        ! Distances to the current abscissa, and their q-th order statistic. Only
        ! that ONE order statistic is needed -- it is the bandwidth -- so it is
        ! selected rather than obtained by sorting all n distances: selection is
        ! O(n) per point where the previous insertion sort was O(n^2) per point on
        ! this input (the distances form a V about x_i, never near-sorted), making
        ! the whole routine cubic in n. dsel is scratch because SelectKth reorders.
        do j = 1, n
            d(j) = abs(x(j) - x(i))
        end do
        dsel = d
        call SelectKth(dsel, n, q)
        dmax = dsel(q)

        ! dmax = 0 means the q-th nearest neighbour coincides with x_i (>= q
        ! duplicated abscissae): the tricube neighbourhood collapses to zero
        ! width, so the local fit is degenerate.
        if (dmax <= 0.0_real64) then
            status = 2
            deallocate(d, dsel)
            return
        end if

        ! Accumulate the weighted normal equations X^T W X beta = X^T W y in the
        ! centred monomial basis [1, dx, dx^2, ...] (dx = x_j - x_i). The
        ! neighbourhood is exactly the tricube support {j : d_j < dmax}; the
        ! q-th neighbour itself lands at d = dmax with weight 0.
        xtwx = 0.0_real64
        xtwy = 0.0_real64
        npos = 0
        do j = 1, n
            if (d(j) < dmax) then
                wj = (1.0_real64 - (d(j) / dmax)**3)**3
                if (wj > 0.0_real64) then
                    npos = npos + 1
                    dx = x(j) - x(i)
                    basis(1) = 1.0_real64
                    do k = 2, p
                        basis(k) = basis(k - 1) * dx
                    end do
                    do k = 1, p
                        do kk = 1, p
                            xtwx(k, kk) = xtwx(k, kk) + wj * basis(k) * basis(kk)
                        end do
                        xtwy(k) = xtwy(k) + wj * basis(k) * y(j)
                    end do
                end if
            end if
        end do

        ! Fewer positive-weight points than free coefficients -> under-determined.
        if (npos < p) then
            status = 2
            deallocate(d, dsel)
            return
        end if

        ! Solve the SPD system; dposv overwrites xtwy with the coefficients.
        ! info > 0 flags a non-positive-definite (rank-deficient) neighbourhood.
        call dposv('U', p, 1, xtwx, p, xtwy, p, info)
        if (info /= 0) then
            status = 2
            deallocate(d, dsel)
            return
        end if

        ! Intercept = fit at the centred origin = fit at x_i.
        ys(i) = xtwy(1)
    end do

    deallocate(d, dsel)
end subroutine LoessFit

end module MTLoess
