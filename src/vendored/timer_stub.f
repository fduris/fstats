c
c  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD
c  License" or "3-clause license"). Please read the attached
c  README.vendored note for the full license text and provenance.
c
c  This file REPLACES the reference distribution's timer.f. The stock
c  timer() reads the process CPU clock (cpu_time); L-BFGS-B uses that
c  only to fill diagnostic timing fields in dsave and never for a
c  convergence decision. Returning a constant 0 makes the optimizer
c  bit-for-bit deterministic across runs and machines, which is what the
c  golden-value regression tests rely on.
c
      subroutine timer(ttime)
      double precision ttime
      ttime = 0.0d0
      return
      end
