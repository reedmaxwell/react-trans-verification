#!/usr/bin/env python
"""
Parameterized ParFlow+Alquimia chemistry runner for the Amanzi benchmark cases.

Builds the common 1D flow/grid/transport deck (identical to the calcite/tracer
cases — 100 m, 100 cells, fixed-head flow, 50 yr) and switches only the
chemistry block per case. Runs ON (install-on) and leaves output under
verification-runs/cases/<case>/out for cross-check against the Amanzi PFLOTRAN
reference.

Usage:
  env PARFLOW_DIR=<install-on> PYTHONPATH=<repo>/pftools/python \
      <subsettools-python> run_case.py <case>
where <case> in {tritium, ion_exchange, surface_complexation, isotherms}.
"""
import os
import sys
import shutil
from parflow import Run
from parflow.tools.fs import cp, mkdir, get_absolute_path

HERE = os.path.dirname(os.path.abspath(__file__))
# Workspace root holding parflow/ and rt-stack/. Override with RT_ROOT;
# the default placeholder must be edited to your local path.
ROOT = os.environ.get("RT_ROOT", "/path/to/parflow_RT_merge")
AMANZI = os.path.join(
    ROOT, "rt-stack/src/amanzi/test_suites/benchmarking/chemistry"
)
JB = os.path.join(ROOT, "rt-stack/src/PF-Alquimia_verification")

# Per-case chemistry: engine deck, files to stage, contaminant slots, CFL,
# and which Chemistry.Print* flags to enable.
CASES = {
    "tritium": dict(
        src=f"{AMANZI}/tritium_1d",
        infile="1d-tritium-crunch.in",
        stage=["1d-tritium-crunch.in", "tritium.dbs", "aqueous.dbs"],
        contaminants="tce dummy",
        cfl=0.6,
        prints=["PrimaryMobile"],
    ),
    "ion_exchange": dict(
        src=f"{JB}/ion_exchange_1d/crunch",
        infile="1d-ion-exchange-crunch.in",
        stage=["1d-ion-exchange-crunch.in", "ion-exchange.dbs"],
        contaminants="tce dummy dummy2 3",
        cfl=0.6,
        prints=["PrimaryMobile", "PrimarySorbed"],
    ),
    "surface_complexation": dict(
        src=f"{JB}/surface_complexation_1d",
        infile="1d-surface-complexation-crunch.in",
        stage=["crunch/1d-surface-complexation-crunch.in", "surface-complexation.dbs"],
        contaminants="tce dummy dummy2 3 4",
        cfl=0.6,
        prints=["PrimaryMobile", "PrimarySorbed", "pH", "SurfSiteDens"],
    ),
    "isotherms": dict(
        src=f"{JB}/isotherms_1d/crunch",
        infile="1d-isotherms-crunch.in",
        stage=["1d-isotherms-crunch.in", "isotherms.dbs"],
        contaminants="tce",
        cfl=0.6,
        prints=["PrimaryMobile", "PrimarySorbed"],
    ),
}


def build_common(r):
    r.FileVersion = 4
    r.Process.Topology.P = 1
    r.Process.Topology.Q = 1
    r.Process.Topology.R = 1

    r.ComputationalGrid.Lower.X = 0.0
    r.ComputationalGrid.Lower.Y = 0.0
    r.ComputationalGrid.Lower.Z = 0.0
    r.ComputationalGrid.DX = 1.0
    r.ComputationalGrid.DY = 1.0
    r.ComputationalGrid.DZ = 1.0
    r.ComputationalGrid.NX = 100
    r.ComputationalGrid.NY = 1
    r.ComputationalGrid.NZ = 1

    r.GeomInput.Names = (
        "domain_input background_input source_region_input concen_region_input"
    )
    r.GeomInput.domain_input.InputType = "Box"
    r.GeomInput.domain_input.GeomName = "domain"
    r.Geom.domain.Lower.X = 0.0
    r.Geom.domain.Lower.Y = 0.0
    r.Geom.domain.Lower.Z = 0.0
    r.Geom.domain.Upper.X = 100.0
    r.Geom.domain.Upper.Y = 1.0
    r.Geom.domain.Upper.Z = 1.0
    r.Geom.domain.Patches = "left right front back bottom top"

    r.GeomInput.background_input.InputType = "Box"
    r.GeomInput.background_input.GeomName = "background"
    r.Geom.background.Lower.X = -99999999.0
    r.Geom.background.Lower.Y = -99999999.0
    r.Geom.background.Lower.Z = -99999999.0
    r.Geom.background.Upper.X = 99999999.0
    r.Geom.background.Upper.Y = 99999999.0
    r.Geom.background.Upper.Z = 99999999.0

    r.GeomInput.source_region_input.InputType = "Box"
    r.GeomInput.source_region_input.GeomName = "source_region"
    r.Geom.source_region.Lower.X = 0.0
    r.Geom.source_region.Lower.Y = 0.0
    r.Geom.source_region.Lower.Z = 0.0
    r.Geom.source_region.Upper.X = 100.0
    r.Geom.source_region.Upper.Y = 1.0
    r.Geom.source_region.Upper.Z = 1.0

    r.GeomInput.concen_region_input.InputType = "Box"
    r.GeomInput.concen_region_input.GeomName = "concen_region"
    r.Geom.concen_region.Lower.X = 0.0
    r.Geom.concen_region.Lower.Y = 0.0
    r.Geom.concen_region.Lower.Z = 0.0
    r.Geom.concen_region.Upper.X = 100.0
    r.Geom.concen_region.Upper.Y = 1.0
    r.Geom.concen_region.Upper.Z = 1.0

    r.Geom.Perm.Names = "background"
    r.Geom.background.Perm.Type = "Constant"
    r.Geom.background.Perm.Value = 0.25
    r.Perm.TensorType = "TensorByGeom"
    r.Geom.Perm.TensorByGeom.Names = "background"
    r.Geom.background.Perm.TensorValX = 1.0
    r.Geom.background.Perm.TensorValY = 1.0
    r.Geom.background.Perm.TensorValZ = 1.0

    r.Phase.Names = "water"
    r.Phase.water.Density.Type = "Constant"
    r.Phase.water.Density.Value = 1.0
    r.Phase.water.Viscosity.Type = "Constant"
    r.Phase.water.Viscosity.Value = 1.0

    r.Gravity = 1.0

    r.TimingInfo.BaseUnit = 1.0
    r.TimingInfo.StartCount = 0
    r.TimingInfo.StartTime = 0.0
    r.TimingInfo.StopTime = 50.0
    r.TimingInfo.DumpInterval = 10.0

    r.Geom.Porosity.GeomNames = "background"
    r.Geom.background.Porosity.Type = "Constant"
    r.Geom.background.Porosity.Value = 0.25

    r.Domain.GeomName = "domain"

    r.Phase.water.Mobility.Type = "Constant"
    r.Phase.water.Mobility.Value = 1.0

    r.Wells.Names = ""

    r.Cycle.Names = "constant"
    r.Cycle.constant.Names = "alltime"
    r.Cycle.constant.alltime.Length = 20
    r.Cycle.constant.Repeat = -1

    r.BCPressure.PatchNames = "left right front back bottom top"
    r.Patch.left.BCPressure.Type = "DirEquilRefPatch"
    r.Patch.left.BCPressure.Cycle = "constant"
    r.Patch.left.BCPressure.RefGeom = "domain"
    r.Patch.left.BCPressure.RefPatch = "bottom"
    r.Patch.left.BCPressure.alltime.Value = 200.0
    r.Patch.right.BCPressure.Type = "DirEquilRefPatch"
    r.Patch.right.BCPressure.Cycle = "constant"
    r.Patch.right.BCPressure.RefGeom = "domain"
    r.Patch.right.BCPressure.RefPatch = "bottom"
    r.Patch.right.BCPressure.alltime.Value = 100.0
    for pn in ("top", "bottom", "back", "front"):
        patch = getattr(r.Patch, pn)
        patch.BCPressure.Type = "FluxConst"
        patch.BCPressure.Cycle = "constant"
        patch.BCPressure.alltime.Value = 0.0

    r.PhaseSources.water.Type = "Constant"
    r.PhaseSources.water.GeomNames = "background"
    r.PhaseSources.water.Geom.background.Value = 0.0

    r.PhaseConcen.water.tce.Type = "Constant"
    r.PhaseConcen.water.tce.GeomNames = "concen_region"
    r.PhaseConcen.water.tce.Geom.concen_region.Value = 0.1

    r.SpecificStorage.Type = "Constant"
    r.SpecificStorage.GeomNames = "background"
    r.Geom.background.SpecificStorage.Value = 1.0e-5

    r.Phase.water.HeatCapacity.Type = "Constant"
    r.Phase.water.HeatCapacity.GeomNames = "background"
    r.Phase.water.Geom.background.HeatCapacity.Value = 4000.0

    for s in ("X", "Y"):
        obj = getattr(r, f"TopoSlopes{s}")
        obj.Type = "Constant"
        obj.GeomNames = "domain"
        obj.Geom.domain.Value = 0.0

    r.Mannings.Type = "Constant"
    r.Mannings.GeomNames = "domain"
    r.Mannings.Geom.domain.Value = 2.3e-7

    # Retardation for the (inert) ParFlow contaminant slots
    r.Geom.Retardation.GeomNames = "background"
    r.Geom.background.tce.Retardation.Type = "Linear"
    r.Geom.background.tce.Retardation.Rate = 0.0


def main():
    case = sys.argv[1]
    p = CASES[case]
    r = Run(f"{case}_pf", __file__)
    # Contaminants first: the PhaseConcen/Retardation children in build_common
    # are created by the Contaminants.Names handler.
    r.Contaminants.Names = p["contaminants"]
    r.Contaminants.tce.Degradation.Value = 0.0
    build_common(r)

    r.Solver.Chemistry = "Alquimia"
    r.Chemistry.Engine = "CrunchFlow"
    r.Chemistry.InputFile = p["infile"]
    r.GeochemCondition.Type = "Constant"
    r.GeochemCondition.GeomNames = "concen_region"
    r.GeochemCondition.Names = "initial"
    r.GeochemCondition.Geom.concen_region.Value = "initial"
    r.BCConcentration.GeochemCondition.Names = "west"
    r.BCConcentration.PatchNames = "left"
    r.Patch.left.BCConcentration.Type = "Constant"
    r.Patch.left.BCConcentration.Value = "west"
    r.Chemistry.ParFlowTimeUnits = "years"
    for flag in p["prints"]:
        setattr(r.Chemistry, f"Print{flag}", True)

    r.Solver.MaxIter = 50000
    r.Solver.CFL = p["cfl"]
    r.Solver.AdvectOrder = 2
    r.Solver.AdvectEnforceMinMax = True
    r.Solver.RelTol = 1.0e-35
    r.Solver.AbsTol = 1.0e-50
    r.Solver.Nonlinear.ResidualTol = 1.0e-15

    out = get_absolute_path(f"{case}/out")
    shutil.rmtree(out, ignore_errors=True)
    mkdir(out)
    for f in p["stage"]:
        cp(f"{p['src']}/{f}", out)
    r.run(working_directory=out)
    print(f"\nDONE {case} -> {out}")


if __name__ == "__main__":
    main()
