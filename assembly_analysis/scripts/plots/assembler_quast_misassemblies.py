#!/usr/bin/env python3
"""Standalone panel d: misassemblies, QUAST-derived.

One of six single-panel scripts (assembler_quast_nga50.py,
_genome_fraction.py, _duplication_ratio.py, _misassemblies.py,
_mismatches.py, _indels.py) that back Figure 4 (de novo assembly
performance) in the thesis report. Drawn through
quast_single_panels.render_panel(), so the grouped-bar style, sample
colors, axis styling and QUAST cross-check are identical across all six
panels -- each bar is one individual GIAB sample observation, never an
average. Saved on its own, with no baked-in panel letter -- see
quast_single_panels.py's docstring for why.

Source: assembly_benchmark_30x.tsv column "misassemblies", cross-checked
against each sample's QUAST report.tsv field "# misassemblies".
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

OUTPUT_PNG = (
    PROJECT / "assembly_analysis" / "figures" / "30x" / "final"
    / "assembler_quast_misassemblies.png"
)


def main() -> int:
    data = render_panel("misassemblies", OUTPUT_PNG, show_strategy=True)
    print(f"\nRows: {len(data)}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PNG.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
