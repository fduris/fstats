! =======================================================================
! MTLoess -- Native loess (locally-weighted polynomial regression) matching R's
! loess(..., family="gaussian") with surface="direct", used to replace R's
! loess() in the plot overlays and the light-mass/speed loess columns.
!
! This is the DOCUMENTED FALLBACK; the primary netlib dloess vendor
! path was declined. dloess' AT&T Bell Labs license IS permissive ("Permission
! to use, copy, modify, and distribute this software for any purpose without
! fee is hereby granted ...", Copyright (c) 1989, 1992 by AT&T), so the license
! gate passes -- but its Fortran core loessf.f FAILS the statelessness gate:
! 29 SAVE statements (per-subroutine execution counters `execnt`, initialised
! once via DATA and incremented every call, never reset; machine constants
! `machin`/`machep` computed only on the first call, `if(execnt.eq.1)`, and
! cached via SAVE for every later call). None are per-call-reinitialised, so
! the spec's stateless-only rule rules the vendor out. Reproducing R's loess()
! would additionally need the C driver stack (loess.c/loessc.c/predict.c: the
! kd-tree interpolation surface), out of scope for this Fortran-only library.
! The fit is therefore done directly here, exactly as R's surface="direct".
!
! Algorithm, per evaluation abscissa x_i (family="gaussian" => NO robustness
! iterations; a single weighted fit):
!   q     = min(n, floor(span*n + Q_TOL))      neighbourhood size
!   dmax  = the q-th smallest |x_j - x_i|       bandwidth (q-th nearest neighbour)
!   w_j   = (1 - (|x_j-x_i|/dmax)^3)^3   for |x_j-x_i| < dmax, else 0   (tricube)
! then fit a weighted polynomial of the requested `degree` (1 or 2) in the
! CENTRED abscissa (x_j - x_i) by weighted least squares; the fitted value at
! x_i is the intercept coefficient (centring makes the intercept the fit AT
! x_i, so no re-evaluation is needed). When q = n the farthest neighbour sits
! at d = dmax and receives weight exactly 0 via (1 - 1^3)^3.
!
! Reproduces R's surface="direct" to ~1e-13 on the goldens (span 0.6 /
! 0.35 / 1.0, degree 1 and 2, and a duplicate-abscissa set). It differs from
! R's DEFAULT surface="interpolate" only by interpolation on R's kd-tree cells;
! that difference is bounded well under 5e-4*sd(y) on the goldens. Stateless:
! no SAVE, no module state, no I/O; all working
! storage is local to the call.
! =======================================================================
module MTLoess

use, intrinsic :: iso_fortran_env, only : real64, int32
use MTSelect, only : SelectKth
implicit none
private

public :: LoessFit

! Rounding tolerance added before truncating span*n to the neighbourhood size,
! matching the `floor(n*span + 1e-5)` of R's loess_workspace (loessc.c). It is
! NOT cosmetic: span*n is generally inexact in binary, so a bare floor drops a
! whole neighbour whenever the product lands a hair below an integer -- e.g.
! span = 0.35, n = 180 gives 62.99999999999999 and hence q = 62 where R uses
! 63. That one neighbour moved the fitted values by 1.6e-3 * sd(y), three times
! the interpolate-vs-direct bound quoted above. Affected (span, n) pairs are
! ordinary ones (span 0.35 with n = 180/340/360, span 0.7 with n = 90/170/180/
! 330..360), so this tolerance is required for R parity, not just tidiness.
real(real64), parameter :: Q_TOL = 1.0e-5_real64

contains

! -----------------------------------------------------------------------
! LoessFit -- Direct-surface loess of `y` on `x` at the input abscissae.
! `degree` is 1
! or 2; `span` in (0, 1] selects the neighbourhood fraction (a span above 1
! is clamped to the full sample by the q = min(n, ...) below, which is NOT
! R's alpha > 1 bandwidth inflation -- the ComputeLoessFit export documents
! (0, 1] as the supported range). On return `ys`
! holds the fitted values and `status` is 0 (computed) or 2 (a local
! neighbourhood was degenerate: fewer than degree+1 positive-weight points,
! all-coincident abscissae in the neighbourhood, or a singular weighted
! normal matrix). R degrades such a neighbourhood to a pseudo-inverse with a
! warning; this routine instead surfaces status=2 rather than emit an
! unreliable smoothed value. Preconditions (n >= 1, degree in {1,2},
! span > 0) are enforced by the ComputeLoessFit export.
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
    ! The Q_TOL before the floor is R's (see the module header).
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
