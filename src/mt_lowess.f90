! =======================================================================
! MTLowess -- Robust locally weighted scatterplot smoothing (LOWESS,
! Cleveland 1979) over a neighbourhood given as a point count, with a
! guard on the robust residual scale.
!
! The fit is fstats_smoothing's lowess with delta = 0, step for step:
!   - each smoothed value is a weighted local line through the ns nearest
!     neighbours, with tricube weights (1 - (r/h)^3)^3 over the distance h
!     to the farthest of them, falling back to the weighted mean when the
!     weighted spread of the neighbourhood is below 0.001 of the x range;
!   - after a pass, cmad = 6 * median|residual| (three times the sum of the
!     two middle order statistics) sets bisquare robustness weights
!     (1 - (r/cmad)^2)^2, with residuals below 0.001 cmad weighted 1 and
!     those above 0.999 cmad weighted 0, for the next pass;
!   - nsteps robustness passes follow the first, plain one.
! Two things differ. The neighbourhood is the count ns itself rather than
! a fraction of n. And the robustness passes stop as soon as cmad falls
! below SCALE_GUARD * mean|y|, keeping the values of the pass just
! computed: once the local fits reproduce the data the residuals are
! rounding noise, and bisquare weights built on noise discard arbitrary
! points and move the curve by amounts that depend on the last digit of
! the input.
!
! Stateless: no SAVE, no module variables, no I/O, all working storage
! local to the call.
! =======================================================================
module MTLowess

use, intrinsic :: iso_fortran_env, only : real64, int32
use linalg, only : sort
implicit none
private

public :: GuardedLowess

! The robust scale below which the robustness passes stop, relative to the
! mean absolute value of the input series y.
real(real64), parameter :: SCALE_GUARD = 1.0e-7_real64

contains

! -----------------------------------------------------------------------
! GuardedLowess -- Robust LOWESS of y on x at the input abscissae.
! x must be non-decreasing. ns is the neighbourhood size in points and is
! clamped into [2, n]; nsteps >= 0 is the largest number of robustness
! passes. A series of fewer than two points is returned unchanged.
! Preconditions (n >= 1, nsteps >= 0) are enforced by the
! ComputeGuardedLowessFit export.
! -----------------------------------------------------------------------
subroutine GuardedLowess(n, x, y, ns, nsteps, ys)
    integer(int32), intent(in) :: n, ns, nsteps
    real(real64), intent(in) :: x(n), y(n)
    real(real64), intent(out) :: ys(n)

    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: three = 3.0d0
    real(real64), parameter :: p001 = 1.0d-3
    real(real64), parameter :: p999 = 0.999d0

    logical :: ok
    integer(int32) :: iter, i, j, nleft, nright, nsize, last, m1, m2
    real(real64) :: d1, d2, denom, alpha, cut, eps, cmad, c1, c9, r, ymean
    real(real64), allocatable :: rw(:), res(:), w(:)

    if (n < 2) then
        ys = y
        return
    end if

    allocate(rw(n), res(n), w(n))

    nsize = max(min(ns, n), 2)
    eps = epsilon(eps)
    ymean = zero
    do i = 1, n
        ymean = ymean + abs(y(i))
    end do
    ymean = ymean / real(n, real64)

    do iter = 1, nsteps + 1
        nleft = 1
        nright = nsize
        last = 0
        i = 1
        do
            do while (nright < n)
                d1 = x(i) - x(nleft)
                d2 = x(nright+1) - x(i)
                if (d1 <= d2) exit
                nleft = nleft + 1
                nright = nright + 1
            end do

            call LocalFit(n, x, y, x(i), ys(i), nleft, nright, w, iter > 1, rw, ok)
            if (.not.ok) ys(i) = y(i)
            if (last < i - 1) then
                denom = x(i) - x(last)
                do j = last + 1, i - 1
                    alpha = (x(j) - x(last)) / denom
                    ys(j) = alpha * ys(i) + (one - alpha) * ys(last)
                end do
            end if
            last = i
            cut = x(last)
            do i = last + 1, n
                if (x(i) > cut) exit
                if (abs(x(i) - x(last)) < eps) then
                    ys(i) = ys(last)
                    last = i
                end if
            end do
            i = max(last + 1, i - 1)

            if (last >= n) exit
        end do

        res = y - ys
        if (iter > nsteps) exit
        rw = abs(res)
        call sort(rw, .true.)
        m1 = 1 + n / 2
        m2 = n - m1 + 1
        cmad = three * (rw(m1) + rw(m2))
        if (cmad < SCALE_GUARD * ymean) exit
        c9 = p999 * cmad
        c1 = p001 * cmad
        do i = 1, n
            r = abs(res(i))
            if (r <= c1) then
                rw(i) = one
            else if (r > c9) then
                rw(i) = zero
            else
                rw(i) = (one - (r / cmad)**2)**2
            end if
        end do
    end do
end subroutine GuardedLowess

! -----------------------------------------------------------------------
! LocalFit -- The weighted local line of one LOWESS pass, evaluated at xs
! over the neighbourhood x(nleft:nright) (fstats_smoothing's lowest).
! With userw the tricube weights are multiplied by the robustness weights
! rw. ok is false when every weight is zero, and ys is then left unset.
! -----------------------------------------------------------------------
subroutine LocalFit(n, x, y, xs, ys, nleft, nright, w, userw, rw, ok)
    integer(int32), intent(in) :: n, nleft, nright
    real(real64), intent(in) :: x(n), y(n), rw(n)
    real(real64), intent(in) :: xs
    real(real64), intent(out) :: ys
    real(real64), intent(inout) :: w(n)
    logical, intent(in) :: userw
    logical, intent(out) :: ok

    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: p001 = 1.0d-3
    real(real64), parameter :: p999 = 0.999d0

    integer(int32) :: j, nrt
    real(real64) :: range, h, h9, h1, a, b, c, r

    range = x(n) - x(1)
    h = max(xs - x(nleft), x(nright) - xs)
    h9 = p999 * h
    h1 = p001 * h
    a = zero

    do j = nleft, n
        w(j) = zero
        r = abs(x(j) - xs)
        if (r <= h9) then
            if (r > h1) then
                w(j) = (one - (r / h)**3)**3
            else
                w(j) = one
            end if
            if (userw) w(j) = rw(j) * w(j)
            a = a + w(j)
        else if (x(j) > xs) then
            exit
        end if
    end do

    nrt = j - 1
    if (a <= zero) then
        ok = .false.
    else
        ok = .true.
        w(nleft:nrt) = w(nleft:nrt) / a
        if (h > zero) then
            a = zero
            do j = nleft, nrt
                a = a + w(j) * x(j)
            end do
            b = xs - a
            c = zero
            do j = nleft, nrt
                c = c + w(j) * (x(j) - a)**2
            end do
            if (sqrt(c) > p001 * range) then
                b = b / c
                do j = nleft, nrt
                    w(j) = w(j) * (one + b * (x(j) - a))
                end do
            end if
        end if
        ys = zero
        do j = nleft, nrt
            ys = ys + w(j) * y(j)
        end do
    end if
end subroutine LocalFit

end module MTLowess
