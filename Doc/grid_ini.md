# PLUTO-style grids in pluto.ini

The `[Grid]` block accepts both the original gPLUTO three-field syntax and the
PLUTO 4.1 patch syntax. No change to the physical evolution or reconstruction
scheme is intended: ini patches supply cell interfaces to the existing grid
and mapped-metric machinery.

```ini
[Grid]
X1-grid 3 -2000 58 s -1024 320 u 1024 58 s 2000
X2-grid 3 -2000 58 s -1024 320 u 1024 58 s 2000
X3-grid 2 3859.2 37 u 4096 442 s 10000
```

Each line is:

```
Xn-grid Npatch x0 N0 type0 x1 N1 type1 ... xNpatch
```

There may be 1–16 contiguous patches. Coordinates must be finite and strictly
increasing; cell counts must be positive integers. Supported patch types:

- `u` / `uniform`: constant physical cell spacing.
- `s` / `stretched`: geometric widths matched to a completed neighbor. The
  solver finds `r > 0` such that `dx*(r + ... + r^N) = patch_length`. It supports
  both increasing and decreasing widths and the uniform `r=1` limit. As in
  PLUTO MakeGrid, the left neighbor is preferred if both are available; patches
  are processed in index order, restarting the scan after each stretched patch.
  A uniform or logarithmic anchor patch is required; an all-stretched axis is
  rejected instead of looping indefinitely.
- `l+`: PLUTO's shifted logarithmic mesh with increasing widths.
- `l-`: the reflected logarithmic mesh with decreasing widths.

For a logarithmic patch `[L,R]`, `a=abs(L)` must be nonzero and
`b=log(1+(R-L)/a)`. At fractional patch index `q`, the edges are
`L+a*expm1(q*b)` for `l+`, and `R-a*expm1((1-q)*b)` for `l-`.
These are PLUTO's formulas, including negative-coordinate patches; they are
not a direct `exp(log(L)+q*log(R/L))` interpolation.

Ghost widths are held equal to the closest physical cell width, matching
PLUTO's grid extension. Inactive axes retain the existing one-cell behavior, including equal endpoints
used by legacy one-dimensional problem files.
`-xres` scales each patch's cell count with the original resolution ratio,
rounding down and retaining at least one cell per patch. Consequently the
resulting total can differ slightly from the requested total for many patches.
Internal grid coarsening requires the factor to divide every patch count;
incompatible counts produce a diagnostic instead of moving a patch interface.

## Definitions and compatibility

Keep `MAPPED_GRID=YES` for nonuniform grids. A user `MappedGrid()` function is
no longer needed for explicit ini patches: the library provides a weak uniform
fallback, overridden by an existing strong user function. An explicit patch
axis always uses its ini specification; a legacy axis still uses the user map
when mapped grids are enabled. Mixing explicit and legacy axes is supported.
The fallback uses the weak-symbol support of GCC, Clang and NVHPC on Linux.
Nonuniform ini grids with `MAPPED_GRID=NO` are rejected because silently using
them with the uniform reconstruction path would be misleading. Uniform patch
grids work with either setting.

Old configuration remains valid:

```ini
X1-grid -2000 436 2000
```

This has exactly three fields and keeps its original mapping behavior. New
configuration has `3*Npatch+2` fields and is recognized by its field count.
Trailing `#` comments are supported. Unknown types, missing fields, zero cell
counts, invalid endpoints and unanchored stretching are diagnosed explicitly.
The parameter reader now allows up to 4095 characters per line (128 lines and
128 tokens per line); this avoids overflowing the former 128-character line
storage when multiple patches are used.

The grid generation path remains host-side and produces the same global arrays
before MPI subdivision. Runtime patch data use a trivially copyable structure
for the existing byte-wise MPI broadcast. MPI and high-order accuracy have not
been validated by the tests below. Piecewise maps can have nonsmooth derivatives
at joins: this change does not establish high-order convergence there.

## Validation

Tests reside in `robust/gbaseline/9_gre/tests/grid_input/`:

- `verify.cpp` embeds the local PLUTO MakeGrid implementation as an independent
  reference. Twelve cases compare all physical cell interfaces, including GRE,
  uniform patches, both logarithmic orientations, negative coordinates,
  chained stretched patches and logarithmic anchors. Additional tests cover
  contracting and constant-width stretching, coarsening and invalid input.
- `integration.py` runs the actual NVHPC/OpenACC GRE executable. It verifies
  legacy/new uniform node and initial-data equivalence, mixed-grid short
  evolution, the production GRE x/y and z axes separately, `-xres` and invalid
  input diagnostics. A separately linked executable with a strong user map
  confirms legacy calls still select that map and explicit ini patches override
  it; results are in `custom_map/results.json`. Separate axes avoid allocating the full 91-million-cell
  production mesh.

Results are recorded in `verify.log` and `integration.json`. Complete production
runs and numerical equivalence of evolved solutions to PLUTO are outside these
grid-generation checks.
