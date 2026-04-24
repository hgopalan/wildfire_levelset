# Regression test: `regrid_finite_check`

## Purpose

Verifies that `regrid_negative_phi` never produces non-finite (`NaN` / `±Inf`)
values in `fine_phi` or `fine_vel` after regridding, even when running in
parallel across multiple MPI ranks.

This test exercises the fail-safes introduced in `src/regrid_negative_phi.H`:

- `amrex::coarsen` used for safe fine→coarse index mapping.
- FAB-box bounds checking before reading coarse arrays.
- `amrex::Math::isfinite` guards replacing bad values with `0.0`.
- Post-transfer `sanitize_multifab` sweep in `regrid_negative_phi`.

## Build

From the repository root:

```bash
git submodule update --init --recursive
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

## Run (serial)

```bash
./build/levelset tests/regrid_finite_check/inputs.i
```

## Run (parallel — recommended to expose race-condition bugs)

```bash
mpirun -n 4 ./build/levelset tests/regrid_finite_check/inputs.i
```

Use any number of MPI ranks ≥ 1 to stress the parallel regridding paths.

## Verification criteria

1. **No crash** — the simulation completes all 20 steps without aborting.
2. **Finite phi_min** — every `phi_min` printed to stdout must be a finite
   number.  Search for `phi_min` in the output and confirm no `-inf` or `nan`
   appears:

   ```bash
   mpirun -n 4 ./build/levelset tests/regrid_finite_check/inputs.i 2>&1 \
       | grep phi_min
   ```

   Expected output (values will differ, but must be finite):

   ```
   Step 0 : phi_min = -3.0 , phi_max = ...
   Step 1 : phi_min = -3.0 , phi_max = ...
   ...
   ```

3. **No non-finite values in fine level** — optionally enable
   `plot_int = 5` and inspect plotfiles with `amrex::VisMF` / `yt` to confirm
   all fine-level cells contain finite values.
