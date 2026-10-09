#!/usr/bin/env python3
"""Standalone panel a: NGA50 (Mb), QUAST-derived.

One of six single-panel scripts (assembler_quast_nga50.py,
_genome_fraction.py, _duplication_ratio.py, _misassemblies.py,
_mismatches.py, _indels.py) that back Figure 4 (de novo assembly
performance) in the thesis report. Drawn through
quast_single_panels.render_panel(), so the grouped-bar style, sample
colors, axis styling and QUAST cross-check are identical across all six
panels -- each bar is one individual GIAB sample observation, never an
average. Saved on its own, with no baked-in panel letter -- see
quast_single_panels.py's docstring for why.

Source: assembly_benchmark_30x.tsv column "nga50" (bp, converted to Mb for
the axis), cross-checked against each sample's QUAST report.tsv field
"NGA50".

NGA50 spans almost three orders of magnitude (GoldRush HiFi ~0.07 Mb to
Flye ONT ~27 Mb), so on a single linear axis everything except Flye ONT
collapses to near-invisible bars. The y-axis is therefore broken: 0-1 Mb
(all GoldRush, Flye HiFi and Verkko values) and 24-28.5 Mb (Flye ONT).
render_panel() refuses to draw if any value would fall inside the break.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "assembly_analysis").is_dir()
)
sys.path.insert(0, str(Path(__file__).resolve().parent))

from quast_single_panels import render_panel

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_quast_nga50.png"


def main() -> int:
    data = render_panel("nga50_mb", OUTPUT_PNG, y_break={
        "lower_max": 1.0, "upper_min": 24.0, "upper_max": 28.5,
        "lower_ticks": [0, 0.25, 0.5, 0.75], "upper_ticks": [26, 28],
    })
    print(f"\nRows: {len(data)}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PNG.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
