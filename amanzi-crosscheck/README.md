# Amanzi cross-check (Phase 5, item 2)

Independent quantitative check of the react_trans port's 1D calcite and tracer
profiles against **Amanzi's shipped PFLOTRAN benchmark reference** for the same
problem decks (the Amanzi chemistry benchmarking suite the GMD 2025 paper draws
from). This upgrades the verification from regression-vs-self (our
`test/chemistry/correct_output` baselines) to comparison against a second,
independent reactive-transport code.

## What is compared
- **Ours**: the port's validated baselines, `test/chemistry/correct_output/*.pfb`
  at dump step 5 (t = 50 yr) — the exact fields the CTest cases verify.
- **Reference**: Amanzi's PFLOTRAN solution of the identical decks, from the
  `.h5` shipped at
  `rt-stack/src/amanzi/test_suites/benchmarking/chemistry/{calcite_1d/pflotran/os,tracer_1d/pflotran}`.
  (Amanzi's calcite/tracer use the SAME `1d-*-crunch.in` + `.dbs` we run.)
- Fields: calcite — Total Ca [mol/L], pH, Calcite volume fraction; tracer —
  Total tracer [mol/L]. **Units identical on both sides; no conversion applied.**

## Result (see `crosscheck_summary.txt`, figures `*_crosscheck.png`)
Front positions agree with PFLOTRAN to well within one grid cell:
- calcite reaction front (Calcite VF): ours 23.02 m vs ref 23.11 m (Δ −0.09 m)
- Ca / pH fronts: Δ +0.49 / +0.43 m (~half a cell)
- tracer front: ours 49.96 m vs ref 49.87 m (Δ +0.09 m)
RMS differences are 1.5–10% of field range, concentrated at the fronts.

The port's fronts are **sharper** than these PFLOTRAN reference runs (clearest
for the tracer and the calcite Ca shoulder near 50 m). This is a scheme
difference, not an error: the port advects with a 2nd-order Godunov scheme +
min/max limiter near Courant 1 (low numerical dispersion), while the PFLOTRAN
reference is more diffusive. The port's sharp tracer front is the near-analytic
Courant-1 result and matches the paper's tracer panel that our baselines already
reproduce. Large pointwise max-relative-Δ values are localized to the steep
fronts / near-zero tails and are expected when comparing two codes.

## Reproduce
```
python compare_amanzi.py
```
(needs a python with h5py + parflow.read_pfb.) Writes the two PNGs and
`crosscheck_summary.txt` here. The Amanzi subtree was pulled with a sparse
checkout of `test_suites/benchmarking/chemistry` only.

## Extending to the other four cases (Phase 5 item 4)
The Amanzi subtree also ships tritium_1d, ion_exchange_1d,
surface_complexation_1d, isotherms_1d with the same `pflotran/*.h5` reference
layout — this script generalizes to them by adding (h5 path, dataset, our pfb)
tuples once those ParFlow cases are built.
