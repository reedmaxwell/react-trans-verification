#!/usr/bin/env python
"""
Three-way Amanzi cross-check at t = 50 yr:
    PFLOTRAN (Amanzi ref, .h5)  vs  CrunchFlow OS3D (Amanzi ref, Tecplot .out)
    vs  this port (PF-Alquimia-CrunchTope, .pfb).

OS3D is the operator-split CrunchFlow variant, so it is the apples-to-apples
geochemistry-engine comparison for our operator-split coupling. Both Amanzi
references ship in the benchmarking/chemistry subtree.

Units: OS3D totcon* are Mol/kgw (== mol/L for these dilute solutions, no
conversion); pH direct; calcite volume* is "Volume %" so /100 -> VF fraction
to match PFLOTRAN Calcite_VF and our MineralVolfx. Mobile species and pH only;
sorbed deferred (separate unit check, see compare_all.py).

Run with a python that has h5py + parflow.read_pfb.
"""
import os
import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from parflow.tools.io import read_pfb

# Workspace root holding parflow/ and rt-stack/. Override with RT_ROOT;
# the default placeholder must be edited to your local path.
ROOT = os.environ.get("RT_ROOT", "/path/to/parflow_RT_merge")
AMANZI = f"{ROOT}/rt-stack/src/amanzi/test_suites/benchmarking/chemistry"
CORRECT = f"{ROOT}/parflow/test/chemistry/correct_output"
CASES_DIR = f"{ROOT}/verification-runs/cases"
OUT = f"{ROOT}/verification-runs/amanzi-crosscheck"
T50 = "Time:  5.00000E+01 y"

# case -> (h5 path, ours dir, run prefix,
#          [(label, our_field, h5_dataset, os3d_relpath, os3d_col, os3d_scale)])
CFG = {
    "tracer": (
        f"{AMANZI}/tracer_1d/pflotran/1d-tracer.h5", CORRECT, "tracer_pf",
        [("Total tracer [mol/L]", "PrimaryMobile.00.tracer", "Total_Tracer [M]",
          "tracer_1d/crunchflow/os3d/totcon5.out", 1, 1.0)],
    ),
    "calcite": (
        f"{AMANZI}/calcite_1d/pflotran/os/1d-calcite.h5", CORRECT, "calcite_pf",
        [("Total Ca [mol/L]", "PrimaryMobile.02.Ca++", "Total_Ca++ [M]",
          "calcite_1d/crunchflow/os3d/totcon5.out", 3, 1.0),
         ("pH", "pH", "pH",
          "calcite_1d/crunchflow/os3d/pH5.out", 1, 1.0),
         ("Calcite vol. fraction", "MineralVolfx.00.Calcite", "Calcite_VF",
          "calcite_1d/crunchflow/os3d/volume5.out", 1, 0.01)],
    ),
    "tritium": (
        f"{AMANZI}/tritium_1d/pflotran/1d-tritium.h5", f"{CASES_DIR}/tritium/out", "tritium_pf",
        [("Total tritium [mol/L]", "PrimaryMobile.00.Tritium", "Total_Tritium [M]",
          "tritium_1d/crunchflow/totcon5.out", 1, 1.0)],
    ),
    "ion_exchange": (
        f"{AMANZI}/ion_exchange_1d/pflotran/1d-ion-exchange.h5",
        f"{CASES_DIR}/ion_exchange/out", "ion_exchange_pf",
        [("Total Na+ [mol/L]", "PrimaryMobile.00.Na+", "Total_Na+ [M]",
          "ion_exchange_1d/crunchflow/totcon5.out", 1, 1.0),
         ("Total Ca++ [mol/L]", "PrimaryMobile.01.Ca++", "Total_Ca++ [M]",
          "ion_exchange_1d/crunchflow/totcon5.out", 2, 1.0),
         ("Total Mg++ [mol/L]", "PrimaryMobile.02.Mg++", "Total_Mg++ [M]",
          "ion_exchange_1d/crunchflow/totcon5.out", 3, 1.0),
         ("Total Cl- [mol/L]", "PrimaryMobile.03.Cl-", "Total_Cl- [M]",
          "ion_exchange_1d/crunchflow/totcon5.out", 4, 1.0)],
    ),
    "surface_complexation": (
        f"{AMANZI}/surface_complexation_1d/pflotran/1d-surface-complexation.h5",
        f"{CASES_DIR}/surface_complexation/out", "surface_complexation_pf",
        [("Total H+ [mol/L]", "PrimaryMobile.00.H+", "Total_H+ [M]",
          "surface_complexation_1d/crunchflow/totcon5.out", 1, 1.0),
         ("Total Na+ [mol/L]", "PrimaryMobile.01.Na+", "Total_Na+ [M]",
          "surface_complexation_1d/crunchflow/totcon5.out", 2, 1.0),
         ("Total NO3- [mol/L]", "PrimaryMobile.02.NO3-", "Total_NO3- [M]",
          "surface_complexation_1d/crunchflow/totcon5.out", 3, 1.0),
         ("Total Zn++ [mol/L]", "PrimaryMobile.03.Zn++", "Total_Zn++ [M]",
          "surface_complexation_1d/crunchflow/totcon5.out", 4, 1.0),
         ("pH", "pH", "pH",
          "surface_complexation_1d/crunchflow/pH5.out", 1, 1.0)],
    ),
    "isotherms": (
        f"{AMANZI}/isotherms_1d/pflotran/1d-isotherms.h5",
        f"{CASES_DIR}/isotherms/out", "isotherms_pf",
        [("Total A [mol/L]", "PrimaryMobile.00.A", "Total_A [M]",
          "isotherms_1d/crunchflow/totcon5.out", 1, 1.0)],
    ),
}


def pflotran_xy(h5_path, ds):
    with h5py.File(h5_path, "r") as f:
        xe = f["Coordinates"]["X [m]"][:]
        xc = 0.5 * (xe[:-1] + xe[1:])
        v = np.squeeze(f[T50][ds][:]).astype(float)
    return xc, v


def ours_xy(ours_dir, prefix, field, dx=1.0):
    v = np.squeeze(read_pfb(f"{ours_dir}/{prefix}.out.{field}.00005.pfb")).astype(float)
    return (np.arange(v.size) + 0.5) * dx, v


def os3d_xy(relpath, col, scale=1.0):
    """Parse a CrunchFlow Tecplot ASCII profile; keep rows whose first token is
    a float (skips '# Units', TITLE, VARIABLES, ZONE headers)."""
    rows = []
    with open(f"{AMANZI}/{relpath}") as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            try:
                x = float(t[0])
                y = float(t[col])
            except (ValueError, IndexError):
                continue
            rows.append((x, y))
    a = np.array(rows, float)
    return a[:, 0], a[:, 1] * scale


def half_level(x, y):
    lo, hi = float(np.min(y)), float(np.max(y))
    if hi - lo < 1e-30:
        return float("nan")
    s = y - 0.5 * (lo + hi)
    idx = np.where(np.diff(np.sign(s)) != 0)[0]
    if idx.size == 0:
        return float("nan")
    i = idx[0]
    return float(x[i] - s[i] * (x[i + 1] - x[i]) / (s[i + 1] - s[i]))


def main():
    lines = ["Amanzi cross-check — primary ref = CrunchFlow OS3D (operator-split, TVD);"
             " PFLOTRAN shown for context (t=50 yr)",
             "Front = half-max location (m); headline = port vs OS3D. "
             "RMS(port-PFLOTRAN) is a like-for-like (mol/L) shape check.", "=" * 78]
    for case, (h5, odir, prefix, fields) in CFG.items():
        lines.append(f"\n{case.upper()}")
        n = len(fields)
        fig, axes = plt.subplots(n, 1, sharex=True, figsize=(8, 2.6 * n + 1))
        if n == 1:
            axes = [axes]
        for ax, (label, ofield, ds, os3d_rel, ocol, oscale) in zip(axes, fields):
            xr, yr = pflotran_xy(h5, ds)
            xo, yo = ours_xy(odir, prefix, ofield)
            xc, yc = os3d_xy(os3d_rel, ocol, oscale)
            yo_i = np.interp(xr, xo, yo)
            rng = float(np.max(yr) - np.min(yr))
            rms = float(np.sqrt(np.mean((yo_i - yr) ** 2)))
            rms_pct = 100 * rms / rng if rng > 1e-30 else float("nan")
            fr, fo, fc = half_level(xr, yr), half_level(xr, yo_i), half_level(xc, yc)
            fronts = (f"front[m]: OS3D {fc:6.2f}  port {fo:6.2f}  (port-OS3D {fo - fc:+.2f})"
                      f" | PFLOTRAN {fr:6.2f}"
                      if np.isfinite(fc) else "front n/a (monotone/flat conservative species)")
            lines.append(f"  {label:24s} {fronts}  [RMS port-PFLO {rms_pct:4.1f}%]")
            # OS3D is the primary reference (operator-split, TVD — matches the
            # paper's choice); the port is the result; PFLOTRAN is de-emphasized
            # context (global-implicit, carries the fix_diffusion numerical
            # dispersion Glenn Hammond/Sergi kept out of the paper figures).
            ax.plot(xc, yc, "-", c="green", lw=3, label="CrunchFlow OS3D (Amanzi ref)")
            ax.plot(xo, yo, "r*", ms=5, label="this port (PF-Alq-Crunch)")
            ax.plot(xr, yr, ":", c="0.5", lw=1.3, label="PFLOTRAN (w/ num. dispersion; context)")
            ax.set_ylabel(label, fontsize=9)
            ax.legend(fontsize=8, loc="best")
        axes[-1].set_xlabel("Distance (m)")
        fig.suptitle(f"{case} — port vs CrunchFlow OS3D (primary ref); PFLOTRAN for context, t=50 yr")
        fig.tight_layout(rect=[0, 0, 1, 0.98])
        fig.savefig(f"{OUT}/{case}_threeway.png", dpi=110)
        plt.close(fig)
    report = "\n".join(lines)
    with open(f"{OUT}/threeway_summary.txt", "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nFigures (*_threeway.png) + threeway_summary.txt in {OUT}")


if __name__ == "__main__":
    main()
