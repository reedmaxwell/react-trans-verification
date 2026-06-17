#!/usr/bin/env python
"""
Six-case Amanzi cross-check (calcite, tracer + the four added cases).

Compares the port's mobile-species profiles at t = 50 yr against Amanzi's
shipped PFLOTRAN reference for the same decks. Mobile species and pH are
directly comparable (mol/L, dimensionless). Sorbed species are NOT compared
here: the PFLOTRAN sorbed datasets are in mol/m^3-bulk while the port's
PrimarySorbed output units are not yet confirmed equal — that needs a separate
unit check (flagged, not assumed).

calcite/tracer "ours" = test/chemistry/correct_output (validated baselines).
the four added cases "ours" = verification-runs/cases/<case>/out (fresh ON runs).

Run with the subsettools env python (h5py + parflow.read_pfb).
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

# case -> (h5 path, ours dir, run prefix, [(label, our_field, h5_dataset), ...])
CFG = {
    "tracer": (
        f"{AMANZI}/tracer_1d/pflotran/1d-tracer.h5", CORRECT, "tracer_pf",
        [("Total tracer [mol/L]", "PrimaryMobile.00.tracer", "Total_Tracer [M]")],
    ),
    "calcite": (
        f"{AMANZI}/calcite_1d/pflotran/os/1d-calcite.h5", CORRECT, "calcite_pf",
        [("Total Ca [mol/L]", "PrimaryMobile.02.Ca++", "Total_Ca++ [M]"),
         ("pH", "pH", "pH"),
         ("Calcite vol. fraction", "MineralVolfx.00.Calcite", "Calcite_VF")],
    ),
    "tritium": (
        f"{AMANZI}/tritium_1d/pflotran/1d-tritium.h5", f"{CASES_DIR}/tritium/out", "tritium_pf",
        [("Total tritium [mol/L]", "PrimaryMobile.00.Tritium", "Total_Tritium [M]")],
    ),
    "ion_exchange": (
        f"{AMANZI}/ion_exchange_1d/pflotran/1d-ion-exchange.h5",
        f"{CASES_DIR}/ion_exchange/out", "ion_exchange_pf",
        [("Total Na+ [mol/L]", "PrimaryMobile.00.Na+", "Total_Na+ [M]"),
         ("Total Ca++ [mol/L]", "PrimaryMobile.01.Ca++", "Total_Ca++ [M]"),
         ("Total Mg++ [mol/L]", "PrimaryMobile.02.Mg++", "Total_Mg++ [M]"),
         ("Total Cl- [mol/L]", "PrimaryMobile.03.Cl-", "Total_Cl- [M]")],
    ),
    "surface_complexation": (
        f"{AMANZI}/surface_complexation_1d/pflotran/1d-surface-complexation.h5",
        f"{CASES_DIR}/surface_complexation/out", "surface_complexation_pf",
        [("Total H+ [mol/L]", "PrimaryMobile.00.H+", "Total_H+ [M]"),
         ("Total Na+ [mol/L]", "PrimaryMobile.01.Na+", "Total_Na+ [M]"),
         ("Total NO3- [mol/L]", "PrimaryMobile.02.NO3-", "Total_NO3- [M]"),
         ("Total Zn++ [mol/L]", "PrimaryMobile.03.Zn++", "Total_Zn++ [M]"),
         ("pH", "pH", "pH")],
    ),
    "isotherms": (
        f"{AMANZI}/isotherms_1d/pflotran/1d-isotherms.h5",
        f"{CASES_DIR}/isotherms/out", "isotherms_pf",
        [("Total A [mol/L]", "PrimaryMobile.00.A", "Total_A [M]")],
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
    lines = ["Amanzi six-case cross-check — port vs PFLOTRAN, t = 50 yr",
             "(mobile species + pH; sorbed deferred pending unit check)", "=" * 64]
    for case, (h5, odir, prefix, fields) in CFG.items():
        lines.append(f"\n{case.upper()}")
        n = len(fields)
        fig, axes = plt.subplots(n, 1, sharex=True, figsize=(8, 2.6 * n + 1))
        if n == 1:
            axes = [axes]
        for ax, (label, ofield, ds) in zip(axes, fields):
            xr, yr = pflotran_xy(h5, ds)
            xo, yo = ours_xy(odir, prefix, ofield)
            yo_i = np.interp(xr, xo, yo)
            rng = float(np.max(yr) - np.min(yr))
            rms = float(np.sqrt(np.mean((yo_i - yr) ** 2)))
            rms_pct = 100 * rms / rng if rng > 1e-30 else float("nan")
            fr, fo = half_level(xr, yr), half_level(xr, yo_i)
            df = fo - fr
            fstr = (f"front ours {fo:6.2f} vs ref {fr:6.2f} m (Δ{df:+.2f})"
                    if np.isfinite(fr) else "front n/a (monotone/flat)")
            lines.append(f"  {label:24s} RMS {rms_pct:5.2f}% of range | {fstr}")
            ax.plot(xr, yr, "-", c="blue", lw=3, label="PFLOTRAN (Amanzi)")
            ax.plot(xo, yo, "r*", ms=5, label="this port")
            ax.set_ylabel(label, fontsize=9)
            ax.legend(fontsize=8, loc="best")
        axes[-1].set_xlabel("Distance (m)")
        fig.suptitle(f"{case} — port vs Amanzi PFLOTRAN, t=50 yr")
        fig.tight_layout(rect=[0, 0, 1, 0.98])
        fig.savefig(f"{OUT}/{case}_crosscheck.png", dpi=110)
        plt.close(fig)
    report = "\n".join(lines)
    with open(f"{OUT}/crosscheck_all_summary.txt", "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nFigures + summary in {OUT}")


if __name__ == "__main__":
    main()
