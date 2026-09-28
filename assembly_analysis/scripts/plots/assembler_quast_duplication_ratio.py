#!/usr/bin/env python3
"""Standalone panel c: duplication ratio, QUAST-derived.

One of six single-panel scripts (assembler_quast_nga50.py,
_genome_fraction.py, _duplication_ratio.py, _misassemblies.py,
_mismatches.py, _indels.py) that back Figure 4 (de novo assembly
performance) in the thesis report. Drawn through
quast_single_panels.render_panel(), so the grouped-bar style, sample
colors, axis styling and QUAST cross-check are identical across all six
panels -- each bar is one individual GIAB sample observation, never an
average. Saved on its own, with no baked-in panel letter -- see
quast_single_panels.py's docstring for why.

Source: assembly_benchmark_30x.tsv column "duplication_ratio",
cross-checked against each sample's QUAST report.tsv field
"Duplication ratio".

Flye HiFi shows an elevated duplication ratio (~1.33) relative to Flye ONT
(~1.01); this figure reports that observed value only -- it does not claim
a confirmed haplotig/purge mechanism, since no further evidence for that
mechanism is in this benchmark. Verkko's ~2.0-2.1 is expected and
qualitatively different: it is a combined ONT+HiFi diploid, haplotype-
resolved assembly, not unpurged duplication (see the dagger note on
Verkko throughout these figures).
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
    / "assembler_quast_duplication_ratio.png"
)


def main() -> int:
    data = render_panel("duplication_ratio", OUTPUT_PNG)
    print(f"\nRows: {len(data)}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PNG.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
