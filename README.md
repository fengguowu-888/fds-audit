# fds-audit

**Static checker for silent geometry failures in [Fire Dynamics Simulator (FDS)](https://pages.nist.gov/fds-smv/) input files.**

FDS *silently* drops, clips, snaps, or neutralizes geometry that falls outside the mesh, sits on solid cells, or is fully occluded — **without any warning in the `.out` file**. The simulation still runs and produces plausible-looking results, so these failures routinely survive peer review and industrial use.

`fds-audit` reconstructs the solid-cell grid from the input and asserts 11 failure modes **before** you run FDS:

```
$ fds-audit model.fds
================================================================================
model.fds
================================================================================
  [i] MESH              demo: 64,000 cells (40×40×40, dx=0.1000×0.1000×0.1000 m)
  [X] E1                OBST cabinet 完全在域外 → FDS 静默丢弃。XB=(4.5, 5.0, 1.0, 3.0, 0.0, 2.0)
```

> Diagnostics are currently emitted in Chinese (the project's working language); an English output pass is planned.

## Why this matters

Three measured facts motivate this tool:

1. **High frequency (C1).** Auditing 925 NIST example inputs, 30.8% contained at least one silent geometry failure; after the tool's own parser fixes, the residual rate is 0.97% — including two genuine bugs in NIST's own examples.
2. **Biased outputs (C2).** In a 40-case study (4 geometries × 9 failure codes × 5 sensors), 6 of 9 failure codes shifted ASET/HRR by ≥5% (12.5%–149.7%) or caused total failure.
3. **Hard to catch by eye (C3).** A mutation benchmark of 77 injected failures across 3 geometries was caught at 100% recall (95% CI [95.3%, 100%]), with 0 false errors on 74 clean cases.

## The 11 checks

| code | level | silent failure caught |
|------|-------|----------------------|
| `E1` containment | ERROR | OBST/VENT fully outside the mesh → silently dropped |
| `E2` clipping | WARN | partially outside → silently clipped |
| `G1` grid conformance | WARN | face not on a grid plane → silently snapped |
| `V1` vent gas path | ERROR | vent with no gas cell on either side → dead opening |
| `V2` vent area | WARN | effective area < 95% of nominal → partial blockage |
| `O1` overlap | WARN | duplicate / fully nested OBST |
| `O2` occlusion | ERROR/WARN | OBST cells all occupied by other solids → inert |
| `D1` devc in solid | ERROR | sensor point inside a solid cell → meaningless reading |
| `B1` boundary audit | WARN | boundary gas cells without an opening (`XB` or `MB` vent) → INERT "fake wall" |
| `R1` resolution | WARN | D\*/dx below threshold (coarse fire resolution) |
| `S1` snapping | INFO | declared → snapped coordinates (discretization report) |

Plus `MESH` (multi-mesh / transformed-grid) and `PARSE` diagnostics.

## Installation

Requires Python ≥ 3.9. No third-party dependencies (standard library only).

```bash
pip install .
# or, for development:
pip install -e .
```

## Usage

### Command line

```bash
fds-audit model.fds            # single file, full report
fds-audit path/to/cases/       # recurse *.fds, summary table
fds-audit a.fds b.fds c.fds    # multiple files, summary table
```

### Python API

```python
from fds_audit import parse_text, check_model

issues = check_model(parse_text(open("model.fds", encoding="utf-8").read()))
for i in issues:
    if i.level == "ERROR":
        print(i.code, i.msg)
```

Validate a model **before writing it to disk** — useful inside parametric generators:

```python
from fds_audit import parse_text, check_model

txt = build_fds_input(params)          # your generator
errs = [i for i in check_model(parse_text(txt)) if i.level == "ERROR"]
if errs:
    raise RuntimeError("silent geometry failure: " + ", ".join(i.code for i in errs))
```

## How it works

`fds-audit` parses the FDS namelist (`&OBST`, `&VENT`, `&HOLE`, `&MESH`, `&DEVC`, `&SURF`, …), reconstructs the solid-cell grid using FDS's own rule (a cell is solid if its *center* lies inside any OBST), then asserts each check against that grid. Uniformly tiled multi-mesh (MPI decomposition) is merged into one logical grid; non-uniform / gapped / overlapping meshes degrade gracefully to a `MESH:WARN` instead of silently guessing. `&VENT MB='XMIN'…` mesh-boundary vents are resolved as full-boundary openings, so the boundary audit does not flag an intentionally open domain boundary.

## Validation

The numbers above come from three studies; scripts and data live alongside the companion paper:

| study | result |
|-------|--------|
| **C3 mutation benchmark** | 77 injected failures, 100% recall (95% CI [95.3%, 100%]); 74 clean cases, 0 false ERROR |
| **C1 NIST audit** | 925 NIST examples: 30.8% → 0.97% ERROR rate after parser fixes; 2 real NIST bugs found |
| **C2 impact** | 40 cases (4×9×5): 6/9 codes bias ASET/HRR ≥5% or cause total failure |

## Limitations (honest)

- Single-mesh and uniformly-tiled multi-mesh only. Transformed grids (`TRNX/Y/Z_ID`) are skipped with a `MESH:WARN` — a Cartesian approximation would be wrong.
- It validates the **input → grid** discretization, not the simulation itself; it does not run FDS.
- Fire resolution (`R1`) uses a coarse ΣHRRPUA·A estimate; treat it as a flag, not a calibrated value.
- Diagnostics are in Chinese for now.

## Citation

If you use `fds-audit` in research, please cite the companion paper (title/DOI to be finalized on acceptance):

> Wu, F., et al. (2026). *Silent geometry failures in Fire Dynamics Simulator: frequency, impact, and an automated checker.* Combustion Science and Technology. (accepted)

A Zenodo DOI for the software will be minted with the first release.

## License

[MIT](LICENSE)
