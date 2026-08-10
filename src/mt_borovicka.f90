! =======================================================================
! MTBorovicka -- Native per-sample Borovicka straight-line radiant
! solver, replacing the straight-line radiant fit.
!
! A meteor's atmospheric trajectory is modelled as a single straight line, a
! point P (Mm) plus a direction R (the radiant, unitless). Each observation is
! a line of sight from a station S (Mm) along a unit direction u. The fit
! minimises the weighted sum of squared perpendicular line-to-line distances
! between the trajectory and every sight (Borovicka 1990, eq. 5).
!
! To remove the two gauge freedoms of the model -- R is defined only up to a
! positive scale (R -> kR) and P only up to a slide along R (P -> P + sR) -- the
! problem is charted on 4 free parameters q. A pivot component k = argmax|R| is
! pinned (R[k] = 1, P[k] = seedPoint[k]) and the remaining freedoms are
!   q = ( P[free1], P[free2], R[free1]/R[k], R[free2]/R[k] ).
! The cost FF4 (the `BorovickaFF4Cost` scalar loop below) and its analytic
! gradient are minimised over q under box bounds by the vendored reference
! L-BFGS-B 3.0 optimizer (`setulb`, reverse communication -- the same driver
! pattern as MTExpFit.f90). A bound-active (pathological) fit re-pivots ONCE
! onto the new argmax and refits; anything still bound-active, non-converged, or
! degenerate returns NaN (the documented non-convergence signal).
!
! Callers must pass maxit = 300: the factr*epsmch stop uses the corrected
! machine epsilon epsilon(one), which needs more iterations for the same
! factr than a looser epsmch would. Do NOT compensate by changing factr.
!
! The cost is a SCALAR per-point accumulation loop, deliberately not a
! vectorized re-association: the converged minimum shifts by ~1e-8 under
! sum reordering.
!
! Stateless: no SAVE, no module variables mutated at call time, no I/O. All
! working storage is local to each call, so a re-pivot in one call cannot leak
! into another, so a re-pivot in one call cannot leak into the next.
! =======================================================================
module MTBorovicka

use, intrinsic :: iso_fortran_env, only : real64, int32
implicit none
private

public :: BorovickaLineDistance, BorovickaFF4Cost, BorovickaFF4Grad
public :: BorovickaMinimize4

! The 4-parameter optimization dimension (P[free1], P[free2], ratios).
integer(int32), parameter :: NQ = 4

! Minimum cross-product magnitude below which two lines are treated as
! parallel. Also the |R| < this -> degenerate radiant test.
real(real64), parameter :: EPS_PARALLEL = 1.0e-30_real64

! Generous bound on the free radiant ratios R[non-k]/R[k].
real(real64), parameter :: RAD_BOUND = 5.0_real64

! Half-width (Mm) of the box on the two free point components, centred on the
! seed point. BorovickaFitOnce builds those bounds from it and OnBound tests the
! fitted point against the same two walls.
real(real64), parameter :: POINT_BOUND = 1.0_real64

! Proximity to a box wall -- +/-RAD_BOUND for a ratio, seed +/-POINT_BOUND for a
! point component -- that counts as "on the bound".
real(real64), parameter :: ONBOUND_TOL = 1.0e-9_real64

! Finite distance sentinel for a parallel sight. Squared by
! the cost to 1e150, large enough to repel the optimiser yet small enough that
! the weighted SUM over all rows stays far below the ~1.8e308 overflow ceiling,
! so the cost never traps on Inf even with many masked rows.
real(real64), parameter :: PARALLEL_SENTINEL = 1.0e75_real64

contains

! -----------------------------------------------------------------------
! FreeIndices -- The two component indices other than the pivot k, in
! increasing order The order fixes which q slot maps to which axis.
! -----------------------------------------------------------------------
pure subroutine FreeIndices(k, f1, f2)
    integer(int32), intent(in) :: k
    integer(int32), intent(out) :: f1, f2

    select case (k)
    case (1)
        f1 = 2; f2 = 3
    case (2)
        f1 = 1; f2 = 3
    case default   ! k == 3
        f1 = 1; f2 = 2
    end select
end subroutine FreeIndices

! -----------------------------------------------------------------------
! SelectPivot -- Index of the largest-magnitude radiant component. maxloc
! returns the FIRST maximum, so ties break to the lowest index.
! -----------------------------------------------------------------------
pure function SelectPivot(radiant) result(k)
    real(real64), intent(in) :: radiant(3)
    integer(int32) :: k

    integer(int32) :: loc(1)

    loc = maxloc(abs(radiant))
    k = loc(1)
end function SelectPivot

! -----------------------------------------------------------------------
! Reconstruct -- Rebuild the point P and radiant R from the 4 free
! parameters q under the pivot chart: the pinned component gets R[k] = 1
! and P[k] = pointFixed, the two free components take q.
! -----------------------------------------------------------------------
pure subroutine Reconstruct(q, k, pointFixed, P, R)
    real(real64), intent(in) :: q(NQ)
    integer(int32), intent(in) :: k
    real(real64), intent(in) :: pointFixed
    real(real64), intent(out) :: P(3), R(3)

    integer(int32) :: f1, f2

    call FreeIndices(k, f1, f2)
    P = 0.0_real64
    R = 0.0_real64
    P(k) = pointFixed
    P(f1) = q(1)
    P(f2) = q(2)
    R(k) = 1.0_real64
    R(f1) = q(3)
    R(f2) = q(4)
end subroutine Reconstruct

! -----------------------------------------------------------------------
! EmCrossU -- g = e_m x u, the cross product of the m-th unit axis with
! the sight direction. Used by the analytic gradient's radiant-ratio
! terms.
! -----------------------------------------------------------------------
pure subroutine EmCrossU(m, ux, uy, uz, g)
    integer(int32), intent(in) :: m
    real(real64), intent(in) :: ux, uy, uz
    real(real64), intent(out) :: g(3)

    select case (m)
    case (1)   ! e1 x u = (0, -u3, u2)
        g(1) = 0.0_real64; g(2) = -uz; g(3) = uy
    case (2)   ! e2 x u = (u3, 0, -u1)
        g(1) = uz; g(2) = 0.0_real64; g(3) = -ux
    case default   ! e3 x u = (-u2, u1, 0)
        g(1) = -uy; g(2) = ux; g(3) = 0.0_real64
    end select
end subroutine EmCrossU

! -----------------------------------------------------------------------
! BorovickaLineDistance -- Signed line-to-line distance between the
! trajectory (point P, direction R) and the sight (station S, unit
! direction u), Borovicka eq. 5. Scale- and translation-gauge invariant:
! unchanged under R -> kR and P -> P + sR. Returns the large finite
! PARALLEL_SENTINEL (never Inf/NaN) when the lines are parallel, so the
! squared cost stays finite.
! -----------------------------------------------------------------------
pure function BorovickaLineDistance(P, R, S, u) result(d)
    real(real64), intent(in) :: P(3), R(3), S(3), u(3)
    real(real64) :: d

    real(real64) :: cr(3), b

    cr(1) = R(2)*u(3) - R(3)*u(2)
    cr(2) = R(3)*u(1) - R(1)*u(3)
    cr(3) = R(1)*u(2) - R(2)*u(1)
    b = sqrt(cr(1)*cr(1) + cr(2)*cr(2) + cr(3)*cr(3))
    if (b < EPS_PARALLEL) then
        d = PARALLEL_SENTINEL
        return
    end if
    d = ((P(1)-S(1))*cr(1) + (P(2)-S(2))*cr(2) + (P(3)-S(3))*cr(3)) / b
end function BorovickaLineDistance

! -----------------------------------------------------------------------
! BorovickaFF4Cost -- FF4 cost: weighted sum of squared line distances
! over all sights. This is the SCALAR per-point accumulation loop,
! deliberately not a vectorized re-association: the converged minimum
! shifts ~1e-8 under sum reordering. A parallel sight contributes
! PARALLEL_SENTINEL^2 = 1e150.
! -----------------------------------------------------------------------
subroutine BorovickaFF4Cost(&
    nStations, stations, nPoints, stIdx, ux, uy, uz, w, q, k, pointFixed, cost &
)
    integer(int32), intent(in) :: nStations, nPoints
    real(real64), intent(in) :: stations(3, nStations)
    integer(int32), intent(in) :: stIdx(nPoints)
    real(real64), intent(in) :: ux(nPoints), uy(nPoints), uz(nPoints), w(nPoints)
    real(real64), intent(in) :: q(NQ)
    integer(int32), intent(in) :: k
    real(real64), intent(in) :: pointFixed
    real(real64), intent(out) :: cost

    real(real64) :: P(3), R(3), uvec(3), d
    integer(int32) :: i, st

    call Reconstruct(q, k, pointFixed, P, R)
    cost = 0.0_real64
    do i = 1, nPoints
        st = stIdx(i)
        uvec(1) = ux(i)
        uvec(2) = uy(i)
        uvec(3) = uz(i)
        d = BorovickaLineDistance(P, R, stations(:, st), uvec)
        cost = cost + w(i) * d * d
    end do
end subroutine BorovickaFF4Cost

! -----------------------------------------------------------------------
! BorovickaFF4Grad -- Closed-form gradient of FF4 w.r.t. the 4 free
! parameters q, accumulated as
! a scalar per-point loop. Shares
! FF4's near-parallel mask: a masked sight gets weight 0 and b = 1, so the
! 1/b, 1/b^3 terms never produce Inf/NaN and the row contributes nothing.
! grad = (d/dP[free1], d/dP[free2], d/dR[free1], d/dR[free2]).
! -----------------------------------------------------------------------
subroutine BorovickaFF4Grad(&
    nStations, stations, nPoints, stIdx, ux, uy, uz, w, q, k, pointFixed, grad &
)
    integer(int32), intent(in) :: nStations, nPoints
    real(real64), intent(in) :: stations(3, nStations)
    integer(int32), intent(in) :: stIdx(nPoints)
    real(real64), intent(in) :: ux(nPoints), uy(nPoints), uz(nPoints), w(nPoints)
    real(real64), intent(in) :: q(NQ)
    integer(int32), intent(in) :: k
    real(real64), intent(in) :: pointFixed
    real(real64), intent(out) :: grad(NQ)

    real(real64) :: P(3), R(3), cr(3), g(3)
    real(real64) :: b, a, b2, b3, D, wi, PSg, crg
    real(real64) :: gp1, gp2, gr1, gr2
    integer(int32) :: i, st, f1, f2

    call Reconstruct(q, k, pointFixed, P, R)
    call FreeIndices(k, f1, f2)

    gp1 = 0.0_real64
    gp2 = 0.0_real64
    gr1 = 0.0_real64
    gr2 = 0.0_real64

    do i = 1, nPoints
        st = stIdx(i)
        cr(1) = R(2)*uz(i) - R(3)*uy(i)
        cr(2) = R(3)*ux(i) - R(1)*uz(i)
        cr(3) = R(1)*uy(i) - R(2)*ux(i)
        b = sqrt(cr(1)*cr(1) + cr(2)*cr(2) + cr(3)*cr(3))
        a = (P(1)-stations(1,st))*cr(1) + (P(2)-stations(2,st))*cr(2) &
            + (P(3)-stations(3,st))*cr(3)
        wi = w(i)
        ! Shared mask with FF4: zero the masked row's weight AND clamp b to 1 so
        ! the 1/b terms stay finite; the whole row then contributes nothing.
        if (b < EPS_PARALLEL) then
            wi = 0.0_real64
            b = 1.0_real64
        end if
        b2 = b*b
        b3 = b2*b
        D = a / b

        ! Point-free gradient: d cost / d P[free] = sum 2 w a cr[free] / b^2.
        gp1 = gp1 + 2.0_real64 * wi * a * cr(f1) / b2
        gp2 = gp2 + 2.0_real64 * wi * a * cr(f2) / b2

        ! Radiant-ratio gradient: sum 2 w D (PSg/b - a crg / b^3), g = e_free x u.
        call EmCrossU(f1, ux(i), uy(i), uz(i), g)
        PSg = (P(1)-stations(1,st))*g(1) + (P(2)-stations(2,st))*g(2) &
            + (P(3)-stations(3,st))*g(3)
        crg = cr(1)*g(1) + cr(2)*g(2) + cr(3)*g(3)
        gr1 = gr1 + 2.0_real64 * wi * D * (PSg / b - a * crg / b3)

        call EmCrossU(f2, ux(i), uy(i), uz(i), g)
        PSg = (P(1)-stations(1,st))*g(1) + (P(2)-stations(2,st))*g(2) &
            + (P(3)-stations(3,st))*g(3)
        crg = cr(1)*g(1) + cr(2)*g(2) + cr(3)*g(3)
        gr2 = gr2 + 2.0_real64 * wi * D * (PSg / b - a * crg / b3)
    end do

    grad(1) = gp1
    grad(2) = gp2
    grad(3) = gr1
    grad(4) = gr2
end subroutine BorovickaFF4Grad

! -----------------------------------------------------------------------
! OnBound -- True if the fit `q` ended pinned on any wall of the box
! BorovickaFitOnce built around the seed `q0`: a free radiant ratio within
! ONBOUND_TOL of +/-RAD_BOUND, or a free point component within
! ONBOUND_TOL of q0 +/-POINT_BOUND. L-BFGS-B reports 'CONVERGENCE' with a
! coordinate sitting on its wall -- the projected gradient of a
! bound-active coordinate clamps to zero -- so a parameter pinned on the
! seed-relative point box is not a trustworthy minimum any more than a
! pinned ratio is. Both are the pathological case that triggers a one-shot
! re-pivot.
! -----------------------------------------------------------------------
pure function OnBound(q, q0) result(res)
    real(real64), intent(in) :: q(NQ), q0(NQ)
    logical :: res

    res = (abs(abs(q(3)) - RAD_BOUND) < ONBOUND_TOL) .or. &
          (abs(abs(q(4)) - RAD_BOUND) < ONBOUND_TOL) .or. &
          (abs(abs(q(1) - q0(1)) - POINT_BOUND) < ONBOUND_TOL) .or. &
          (abs(abs(q(2) - q0(2)) - POINT_BOUND) < ONBOUND_TOL)
end function OnBound

! -----------------------------------------------------------------------
! BorovickaFitOnce -- One L-BFGS-B minimisation of FF4 under box bounds.
! Bounds: the point-free parameters are q0(1:2) +/-POINT_BOUND Mm, the
! radiant ratios +/-RAD_BOUND. The reverse-communication driver and
! convergence-code mapping match MTExpFit:
!   0 converged, 1 iteration limit, 51 warning, 52 abnormal/other.
! status is 0 when the optimiser ran, 1 when it could not be started at all.
! Stateless: all L-BFGS-B working storage is local to this call.
! -----------------------------------------------------------------------
subroutine BorovickaFitOnce(&
    nStations, stations, nPoints, stIdx, ux, uy, uz, w, k, pointFixed, q0, maxit, lmm, factr, &
    pgtol, qout, convergence, status &
)
    integer(int32), intent(in) :: nStations, nPoints
    real(real64), intent(in) :: stations(3, nStations)
    integer(int32), intent(in) :: stIdx(nPoints)
    real(real64), intent(in) :: ux(nPoints), uy(nPoints), uz(nPoints), w(nPoints)
    integer(int32), intent(in) :: k
    real(real64), intent(in) :: pointFixed
    real(real64), intent(in) :: q0(NQ)
    integer(int32), intent(in) :: maxit, lmm
    real(real64), intent(in) :: factr, pgtol
    real(real64), intent(out) :: qout(NQ)
    integer(int32), intent(out) :: convergence, status

    ! setulb reverse-communication state (written by setulb, opaque here)
    character(len=60) :: task, csave
    logical :: lsave(4)
    integer(int32) :: isave(44)
    real(real64) :: dsave(29)

    ! problem description and working storage for setulb
    integer(int32) :: nbd(NQ)
    integer(int32) :: iwa(3*NQ)
    real(real64), allocatable :: wa(:)
    real(real64) :: lower(NQ), upper(NQ), g(NQ), f
    integer(int32) :: iter, ierr

    ! iprint < 0 suppresses the L-BFGS-B iteration report and the iterate.dat
    ! file. It does NOT cover two unguarded write(6,*) calls in the reference
    ! sources; those are commented out in the vendored copy (see README.vendored),
    ! without which this would print to stdout from inside the library.
    integer(int32), parameter :: IPRINT = -1

    convergence = 52
    status = 0
    qout = q0
    f = 0.0_real64
    g = 0.0_real64

    if (nPoints < 1 .or. nStations < 1 .or. lmm < 1 .or. maxit < 0) then
        status = 1
        return
    end if

    ! setulb workspace: wa(2*m*n + 5*n + 11*m*m + 8*m), iwa(3*n), n=NQ, m=lmm.
    allocate(wa(2*lmm*NQ + 5*NQ + 11*lmm*lmm + 8*lmm), stat=ierr)
    if (ierr /= 0) then
        status = 1
        return
    end if

    ! Point-free parameters are bounded +/-POINT_BOUND Mm about the seed, ratios
    ! +/-RAD_BOUND.
    lower(1) = q0(1) - POINT_BOUND
    upper(1) = q0(1) + POINT_BOUND
    lower(2) = q0(2) - POINT_BOUND
    upper(2) = q0(2) + POINT_BOUND
    lower(3) = -RAD_BOUND
    upper(3) = RAD_BOUND
    lower(4) = -RAD_BOUND
    upper(4) = RAD_BOUND

    ! nbd(i) = 2: every parameter has both a lower and an upper bound.
    nbd = 2
    iter = 0
    task = 'START'

    ! Reverse-communication loop.
    do
        call setulb(&
            NQ, lmm, qout, lower, upper, nbd, f, g, factr, pgtol, wa, iwa, task, IPRINT, csave, &
            lsave, isave, dsave &
        )

        if (task(1:2) == 'FG') then
            call BorovickaFF4Cost(&
                nStations, stations, nPoints, stIdx, ux, uy, uz, w, qout, k, pointFixed, f &
            )
            call BorovickaFF4Grad(&
                nStations, stations, nPoints, stIdx, ux, uy, uz, w, qout, k, pointFixed, g &
            )
        else if (task(1:5) == 'NEW_X') then
            ! Callers must pass maxit = 300 (see the module header).
            iter = iter + 1
            if (iter > maxit) then
                convergence = 1
                exit
            end if
        else if (task(1:4) == 'WARN') then
            convergence = 51
            exit
        else if (task(1:4) == 'CONV') then
            convergence = 0
            exit
        else
            ! 'ERROR...', 'ABNORMAL...', or any unexpected terminal task.
            convergence = 52
            exit
        end if
    end do

    deallocate(wa)
end subroutine BorovickaFitOnce

! -----------------------------------------------------------------------
! BorovickaMinimize4 -- The full per-sample solver. Charts on
! the seed radiant's argmax pivot, fits, and -- if a bound is active -- re-
! pivots once onto the fitted radiant's argmax and refits. Any non-
! convergence, a still-bound-active or same-pivot re-pivot, an all-zero or
! non-finite SEED radiant, or a degenerate (|R| < EPS_PARALLEL) or non-finite
! FITTED radiant yields NaN outputs with converged = 0 (the
! documented non-convergence signal). On success outRadiant is the unit
! radiant and outPoint is the point in Mm (the Pascal wrapper scales to m).
!
! status: 0 the solver ran (converged is 0 or 1); 1 it could not be started
! (bad sizes, an out-of-range station index, or an L-BFGS-B workspace
! allocation failure).
! -----------------------------------------------------------------------
subroutine BorovickaMinimize4(&
    nStations, stations, nPoints, stIdx, ux, uy, uz, w, seedPoint, seedRadiant, maxit, lmm, &
    factr, pgtol, outRadiant, outPoint, converged, status &
)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_is_finite
    integer(int32), intent(in) :: nStations, nPoints
    real(real64), intent(in) :: stations(3, nStations)
    integer(int32), intent(in) :: stIdx(nPoints)
    real(real64), intent(in) :: ux(nPoints), uy(nPoints), uz(nPoints), w(nPoints)
    real(real64), intent(in) :: seedPoint(3), seedRadiant(3)
    integer(int32), intent(in) :: maxit, lmm
    real(real64), intent(in) :: factr, pgtol
    real(real64), intent(out) :: outRadiant(3), outPoint(3)
    integer(int32), intent(out) :: converged, status

    real(real64) :: q0(NQ), qfit(NQ), P(3), R(3)
    real(real64) :: pointFixed, pf2, nrm, nan
    integer(int32) :: k, k2, f1, f2, conv, stat

    nan = ieee_value(0.0_real64, ieee_quiet_nan)
    outRadiant = nan
    outPoint = nan
    converged = 0
    status = 0

    if (nPoints < 1 .or. nStations < 1 .or. lmm < 1 .or. maxit < 0) then
        status = 1
        return
    end if
    ! Every sight must reference a station that exists.
    if (any(stIdx < 1) .or. any(stIdx > nStations)) then
        status = 1
        return
    end if

    ! A non-finite seed radiant cannot be charted (Inf/Inf and NaN/anything both
    ! give NaN ratios, and maxloc's choice of pivot is not even meaningful), so
    ! reject it before the pivot is chosen.
    if (.not. all(ieee_is_finite(seedRadiant))) return   ! NaN outputs, converged = 0

    ! Seed chart: pivot on the largest seed-radiant component.
    k = SelectPivot(seedRadiant)
    pointFixed = seedPoint(k)

    ! The pivot component is the divisor of both ratios below and is the largest
    ! in magnitude, so it is zero only when the entire seed radiant is zero: that
    ! is the one seed that would form 0/0 and hand the optimizer an NaN start.
    ! Reject it here rather than relying on setulb to bounce the NaN back into the
    ! non-convergence path (which it does today, but by its own NaN handling
    ! rather than by contract).
    !
    ! Deliberately an EXACT-ZERO test, not a magnitude threshold: the chart
    ! divides by its own pivot and is therefore scale-free, so an arbitrarily
    ! small nonzero seed is perfectly well determined. (1e-40, 0, 2e-40) has the
    ! same ratios as (0.5, 0, 1) and must still fit -- an |R| < EPS_PARALLEL gate
    ! here would reject a valid seed.
    if (seedRadiant(k) == 0.0_real64) return
    call FreeIndices(k, f1, f2)
    q0(1) = seedPoint(f1)
    q0(2) = seedPoint(f2)
    q0(3) = seedRadiant(f1) / seedRadiant(k)
    q0(4) = seedRadiant(f2) / seedRadiant(k)

    call BorovickaFitOnce(&
        nStations, stations, nPoints, stIdx, ux, uy, uz, w, k, pointFixed, q0, maxit, lmm, factr, &
        pgtol, qfit, conv, stat &
    )
    if (stat /= 0) then
        status = 1
        return
    end if
    if (conv /= 0) return   ! NaN outputs, converged = 0

    call Reconstruct(qfit, k, pointFixed, P, R)

    if (OnBound(qfit, q0)) then
        ! Bound-active fit: re-pivot ONCE onto the reconstructed radiant's argmax
        ! and refit, using purely local k2/pf2. Note pf2
        ! and q02 are built from the RECONSTRUCTED P and R, and R is normalised by
        ! its own k2 component.
        k2 = SelectPivot(R)
        if (k2 == k) return
        call FreeIndices(k2, f1, f2)
        pf2 = P(k2)
        q0(1) = P(f1)
        q0(2) = P(f2)
        q0(3) = R(f1) / R(k2)
        q0(4) = R(f2) / R(k2)

        call BorovickaFitOnce(&
            nStations, stations, nPoints, stIdx, ux, uy, uz, w, k2, pf2, q0, maxit, lmm, factr, &
            pgtol, qfit, conv, stat &
        )
        if (stat /= 0) then
            status = 1
            return
        end if
        if (conv /= 0) return
        call Reconstruct(qfit, k2, pf2, P, R)
        ! q0 now holds the re-pivot seed, so this tests the second fit's own box.
        if (OnBound(qfit, q0)) return   ! still bound-active -> non-converged
    end if

    ! Degenerate-radiant gate, written as .not. (nrm >= EPS) rather than
    ! (nrm < EPS) so an NaN is REJECTED. A NaN fails every ordered comparison, so
    ! the naive form let one fall through and be reported as converged = 1 with
    ! NaN outputs -- the exact opposite of this routine's contract, under which
    ! NaN means non-convergence.
    nrm = sqrt(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))
    if (.not. (nrm >= EPS_PARALLEL)) return   ! degenerate or non-finite radiant

    outRadiant(1) = R(1) / nrm
    outRadiant(2) = R(2) / nrm
    outRadiant(3) = R(3) / nrm
    outPoint(1) = P(1)
    outPoint(2) = P(2)
    outPoint(3) = P(3)
    converged = 1
end subroutine BorovickaMinimize4

end module MTBorovicka
