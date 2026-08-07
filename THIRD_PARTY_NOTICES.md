# Third-party notices

`libfstats` embeds the third-party code listed here.
The 3-clause BSD license of that code requires, for **binary** redistribution, that the copyright notice, the list of conditions and the disclaimer be reproduced "in the documentation and/or other materials provided with the distribution".

**Anything that ships a built `libfstats.so` / `libfstats.dll` must ship this file (or its contents) alongside the binary.**
That includes the compiled libraries checked in at the root of the consuming MT repository, and any release archive or installer built from them.
The notice travelling only inside `src/vendored/` is not sufficient once the binary is distributed on its own.

## L-BFGS-B 3.0

Bound-constrained limited-memory quasi-Newton optimizer, used by the exp-model fit and the Borovicka radiant solver.

- Vendored at `src/vendored/lbfgsb.f`, with the LINPACK subset it needs at `src/vendored/linpack_subset.f`.
- Authors: Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis Morales.
- License: 3-clause ("New"/"Modified") BSD.
  Full text, together with the attribution note explaining that upstream ships the template with unfilled placeholders, is in **`src/vendored/License.txt`**.
- Provenance, the per-file change list, and the deviations from the reference sources are recorded in `src/vendored/README.vendored`.

`src/vendored/timer_stub.f` replaces upstream's `timer.f` and is original to this repository; it carries the same notice by association.

## Not embedded

BLAS and LAPACK are linked, not vendored, and carry their own notices from whichever implementation the build resolves them against.
