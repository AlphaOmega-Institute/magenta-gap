#!/usr/bin/env python3
"""Run every independent check and store its printed output under outputs/.

Usage    python research/independent-checks/run_all.py            (from the repository root, or anywhere)
Needs    Python 3.11+, NumPy 2.x, SciPy 1.10+.  About five minutes in total.
Reports  One line per script: exit status and run time.  Exit status 1 if any script fails to run.

These scripts print the quantities quoted in the November 2026 audit; most do not assert.  A clean run shows that each
script executes and lets the stored outputs be compared by eye or by diff.  Lines that report run times differ between runs.
"""
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
RUNS = [
    ("hull_and_purple_edge.py", [], None),
    ("full_block_fisher.py", [], None),
    ("chart_comparison.py", [], None),
    ("symmetric_body_and_order_interval.py", [], None),
    ("linear_channel_symmetry.py", [], None),
    ("span_collinearity_hilbert.py", [], None),
    ("bounded_reflectance_body.py", [], None),
    ("widest_point_types.py", [], None),
    ("crease_slopes_and_near_black.py", [], None),
    ("continuum_rim.py", [], "continuum_rim_390_700"),
    ("continuum_rim.py", ["830"], "continuum_rim_390_830"),
    ("three_transition_points.py", [], None),
    ("crease_fan_and_550nm.py", [], None),
    ("cubic_law_and_september_numbers.py", [], None),
    ("certificate_review.py", [], None),
    ("flat_face_counterexample.py", [], None),
    ("timelike_direction_step.py", [], None),
    ("effect_body_and_vector_identity.py", [], None),
    ("hue_projectors_holonomy_boost.py", [], None),
    ("white_point_and_chord_region.py", [], None),
    ("ground_plane_interval.py", [], None),
    ("lg_effect_family.py", [], None),
    ("two_encoding_routes.py", [], None),
    ("defect_extensions_and_altered_points.py", [], None),
]


def main():
    out_dir = HERE / "outputs"
    out_dir.mkdir(exist_ok=True)
    failed = 0
    for script, args, name in RUNS:
        start = time.time()
        proc = subprocess.run([sys.executable, str(HERE / script), *args], capture_output=True, text=True)
        text = proc.stdout + (("\n[stderr]\n" + proc.stderr) if proc.stderr.strip() else "")
        (out_dir / f"{name or script[:-3]}.txt").write_text(text)
        status = "ok  " if proc.returncode == 0 else "FAIL"
        failed += proc.returncode != 0
        print(f"{status} {script} {' '.join(args):4s} {time.time() - start:6.1f}s")
    print(f"\n{len(RUNS) - failed} of {len(RUNS)} ran cleanly")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
