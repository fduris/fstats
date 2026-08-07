! =======================================================================
! MTSelect -- Order-statistic primitives shared by the MT numerics: a
! deterministic
! selection routine and the R-compatible median / MAD built on it.
!
! These live in their own module because two unrelated consumers need the same
! primitive -- MTRobust for the medians the mblm estimators are defined by, and
! MTLoess for the q-th nearest-neighbour distance that sets its bandwidth --
! and neither is a natural owner for the other.
!
! Stateless and pure throughout: no SAVE, no module variables, no I/O, all
! working storage in caller-passed arrays or locals.
! =======================================================================
module MTSelect

use, intrinsic :: iso_fortran_env, only : real64, int32
implicit none
private

public :: SelectKth, MedianInplace, MadInplace, MedianOf, MadOf

contains

! -----------------------------------------------------------------------
! SelectKth -- Partially reorders a(1:m) so that a(k) becomes the k-th
! smallest element,
! every earlier element is <= it and every later element is >= it.
!
! Quickselect with a median-of-three pivot: O(m) expected, iterative (no
! recursion), no working storage, and fully DETERMINISTIC -- no
! randomisation, so equal input always yields a bitwise equal result, which
! the reentrancy suite's bitwise-equality requirement depends on.
! -----------------------------------------------------------------------
pure subroutine SelectKth(a, m, k)
    integer(int32), intent(in) :: m, k
    real(real64), intent(inout) :: a(m)

    integer(int32) :: lo, hi, i, j, mid
    real(real64) :: pivot, t

    lo = 1
    hi = m
    do while (lo < hi)
        ! Median-of-three (lo, mid, hi), which also leaves a(lo) <= pivot <= a(hi)
        ! and so bounds both partition scans below without an explicit index test.
        mid = lo + (hi - lo) / 2
        if (a(mid) < a(lo)) then
            t = a(mid); a(mid) = a(lo); a(lo) = t
        end if
        if (a(hi) < a(lo)) then
            t = a(hi); a(hi) = a(lo); a(lo) = t
        end if
        if (a(hi) < a(mid)) then
            t = a(hi); a(hi) = a(mid); a(mid) = t
        end if
        pivot = a(mid)

        ! Hoare partition about the pivot VALUE: on exit a(lo:j) <= pivot and
        ! a(j+1:hi) >= pivot. Both indices always advance across a swap, so the
        ! inner loop terminates even when every element equals the pivot.
        i = lo
        j = hi
        do
            do while (a(i) < pivot)
                i = i + 1
            end do
            do while (a(j) > pivot)
                j = j - 1
            end do
            if (i >= j) exit
            t = a(i); a(i) = a(j); a(j) = t
            i = i + 1
            j = j - 1
        end do

        ! Keep only the side that still contains position k; k stays in [lo, hi],
        ! so the loop ends with lo = hi = k.
        if (k <= j) then
            hi = j
        else
            lo = j + 1
        end if
    end do
end subroutine SelectKth

! -----------------------------------------------------------------------
! MedianInplace -- R's median of a(1:m) -- the central order statistic
! for odd m, the mean of
! the two central ones for even m -- computed by selection. REORDERS `a`, so
! the caller must be done with its order. Identical value to sorting: the
! order statistics of a multiset do not depend on how they are found.
! (A subroutine, not a function: a pure function may not take intent(inout).)
! -----------------------------------------------------------------------
pure subroutine MedianInplace(a, m, med)
    integer(int32), intent(in) :: m
    real(real64), intent(inout) :: a(m)
    real(real64), intent(out) :: med

    integer(int32) :: h

    if (mod(m, 2) == 1) then
        h = (m + 1) / 2
        call SelectKth(a, m, h)
        med = a(h)
    else
        h = m / 2
        call SelectKth(a, m, h)
        ! SelectKth leaves every element past h >= a(h), so the (h+1)-th order
        ! statistic is simply the smallest of that tail.
        med = 0.5_real64 * (a(h) + minval(a(h+1:m)))
    end if
end subroutine MedianInplace

! -----------------------------------------------------------------------
! MadInplace -- R's mad(a) = 1.4826 * median(|a - center|) with center =
! median(a), which
! every caller here has already computed. OVERWRITES `a` with the absolute
! deviations, so one buffer serves both the median and the MAD.
! -----------------------------------------------------------------------
pure subroutine MadInplace(a, m, center, res)
    integer(int32), intent(in) :: m
    real(real64), intent(inout) :: a(m)
    real(real64), intent(in) :: center
    real(real64), intent(out) :: res

    real(real64) :: med

    a = abs(a - center)
    call MedianInplace(a, m, med)
    res = 1.4826_real64 * med
end subroutine MadInplace

! -----------------------------------------------------------------------
! MedianOf -- Non-destructive R median, for the O(n)-sized vectors where
! a private copy
! costs nothing. The big pairwise buffer uses MedianInplace instead.
! -----------------------------------------------------------------------
pure function MedianOf(a, m) result(med)
    integer(int32), intent(in) :: m
    real(real64), intent(in) :: a(m)
    real(real64) :: med

    real(real64) :: tmp(m)

    tmp = a
    call MedianInplace(tmp, m, med)
end function MedianOf

! -----------------------------------------------------------------------
! MadOf -- Non-destructive R mad(a) = 1.4826 * median(|a - median(a)|), for the
! O(n)-sized vectors.
! -----------------------------------------------------------------------
pure function MadOf(a, m) result(res)
    integer(int32), intent(in) :: m
    real(real64), intent(in) :: a(m)
    real(real64) :: res

    real(real64) :: tmp(m)

    tmp = a
    call MadInplace(tmp, m, MedianOf(a, m), res)
end function MadOf

end module MTSelect
