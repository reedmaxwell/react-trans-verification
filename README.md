# react_trans verification & benchmarking — how-to

This directory holds the **verification artifacts** for the ParFlow react_trans
(Alquimia/CrunchFlow reactive-transport) port: the scripts that run the 1D
benchmark cases, the cross-check against Amanzi's PFLOTRAN reference, the output
figures/summaries, and the team report.

It is **analysis, not part of the ParFlow repo** — it lives in a project
workspace one level above the `parflow/` git checkout and is tracked as its own
standalone repo (`reedmaxwell/react-trans-verification`). See
`benchmark-report/react_trans_benchmark_report.pdf` for the writeup; this file is
the operational how-to to reproduce it.

> **Paths:** the scripts here reference a workspace that holds `parflow/`,
> `rt-stack/`, and this tree side by side. Set `RT_ROOT` to that workspace (or
> edit the `/path/to/parflow_RT_merge` placeholder in each script) before
> running — no machine-specific absolute paths are committed.

---

## What's here (inside `verification-runs/`)

```
cases/
  run_case.py              parameterized ParFlow+Alquimia runner for the 4 added cases
  <case>/out/              fresh ON run output per case (tritium, ion_exchange,
                           surface_complexation, isotherms)
amanzi-crosscheck/
  compare_all.py           six-case cross-check vs PFLOTRAN -> figures + summary
  compare_amanzi.py        original calcite+tracer detailed cross-check
  crosscheck_all_summary.txt   numeric results (front positions, RMS) — all six
  *_crosscheck.png         per-case overlay figures
  README.md                cross-check method notes
benchmark-report/
  react_trans_benchmark_report.{md,pdf,html}   the team report (Sergi, Steve)
  figures/                 figures embedded in the report
fig2-calcite/ fig2-tracer/ td-calcite/   earlier acceptance-gate run artifacts
```

## What it depends on (OUTSIDE `verification-runs/`)

| Dependency | Location | What it is |
|---|---|---|
| ParFlow ON build | `../install-on` | ParFlow built with `PARFLOW_ENABLE_ALQUIMIA=ON` |
| Engine stack | `../rt-stack/install/` | OpenMPI 4.1.6, PETSc, CrunchTope, alquimia (+ gfortran-rt) |
| Amanzi PFLOTRAN refs | `../rt-stack/src/amanzi/test_suites/benchmarking/chemistry` | sparse checkout; the independent reference `.h5` |
| jbeisman drivers | `../rt-stack/src/PF-Alquimia_verification` | ParFlow `.tcl` drivers + decks for the 4 added cases |
| Committed baselines | `../parflow/test/chemistry/correct_output` | calcite/tracer "ours" side (validated) |
| Python env | `subsettools` conda env | has `h5py`, `numpy`, `matplotlib`, `parflow.read_pfb` |

Paths below assume `ROOT=/path/to/parflow_RT_merge` (the workspace dir holding
both `parflow/` and this `verification-runs/`). The engine
stack is prebuilt; rebuild recipes (if needed) are
`../docs/build_alquimia_crunch_macos.sh` and
`../docs/rebuild_deps_openmpi416_macos.sh`.

---

## Prerequisites

```sh
ROOT=/path/to/parflow_RT_merge
PF=$ROOT/parflow                       # the git checkout
PYBIN=python   # any python with h5py + numpy + matplotlib + parflow.read_pfb
```

1. A ParFlow ON build at `$ROOT/install-on` (Alquimia + CrunchTope engine).
2. The `subsettools` python (h5py + parflow.read_pfb).
3. The Amanzi reference subtree (one-time sparse checkout):

```sh
cd $ROOT/rt-stack/src
git clone --no-checkout --depth 1 --filter=blob:none \
    https://github.com/amanzi/amanzi.git amanzi
cd amanzi
git sparse-checkout init --cone
git sparse-checkout set test_suites/benchmarking/chemistry
git checkout
```

---

## Replicate the runs

The four added cases (calcite + tracer are already committed regression tests in
the repo; their "ours" baselines live in `test/chemistry/correct_output`):

```sh
cd $ROOT/verification-runs/cases
for c in tritium ion_exchange surface_complexation isotherms; do
  env PARFLOW_DIR=$ROOT/install-on \
      PYTHONPATH=$PF/pftools/python \
      $PYBIN run_case.py $c
done
# output lands in cases/<case>/out/<case>_pf.out.*.pfb
```

`run_case.py` builds the common 1D flow/grid deck (identical to calcite/tracer)
and switches only the per-case chemistry (engine deck, species, print flags) via
its `CASES` dict. It stages the engine `.in`/`.dbs` into the working dir and runs
ParFlow ON.

## Run the cross-check vs PFLOTRAN

```sh
cd $ROOT/verification-runs/amanzi-crosscheck
$PYBIN compare_all.py        # six cases -> *_crosscheck.png + crosscheck_all_summary.txt
$PYBIN compare_amanzi.py     # calcite+tracer, more detailed metrics
```

`compare_all.py` reads our pfb output (calcite/tracer from the committed
`correct_output`; the four added cases from `cases/<case>/out`), extracts the
matching field from the Amanzi PFLOTRAN `.h5` at t=50 yr, and reports front
position + RMS per mobile species/pH. **Units are mol/L (mobile), dimensionless
(pH), mol/m^3 (sorbed) on both sides — no conversion is applied; sorbed units
were verified equal, not assumed.**

## Regenerate the report

```sh
cd $ROOT/verification-runs/benchmark-report
pandoc react_trans_benchmark_report.md -s --self-contained \
    --metadata title="ParFlow react_trans — 1D Benchmark Verification" \
    -o react_trans_benchmark_report.html
pandoc react_trans_benchmark_report.md --pdf-engine=xelatex \
    -V geometry:margin=1in -o react_trans_benchmark_report.pdf
```

---

## Notes / gotchas

- **`PARFLOW_DIR` must point at the ON install** (`install-on`), not the shell
  profile default. The runner launches ParFlow from there.
- **`PYTHONPATH=$PF/pftools/python`** is only needed so the runner picks up the
  in-tree pftools (with the chemistry keys). For the cross-check, `subsettools`'
  own `parflow.read_pfb` is sufficient.
- **Species notation in the report is ASCII** (Ca++, Cl-, NO3-) so the xelatex
  PDF font renders cleanly; superscripts / `≈` / `★` / `m³` break the default
  LaTeX font.
- The cross-check compares **mobile species + pH** quantitatively; sorbed
  species are spot-verified (units + magnitude) but not plotted per-species.
- This whole tree is **uncommitted analysis**. Committing the four cases as
  CTest regression tests (Tier 1) and landing a standalone benchmark suite
  (Tier 2) awaits the team's placement/scope decision — see report §7 and §9.
