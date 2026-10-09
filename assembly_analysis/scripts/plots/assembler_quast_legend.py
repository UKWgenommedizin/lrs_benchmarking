#!/usr/bin/env python3
"""Shared sample legend (HG002/HG003/HG004) for the six QUAST panels of
Figure 4. Drawn once as its own strip and placed above the panels in LaTeX,
mirroring fig3_legend_30x for Figure 3 -- see quast_single_panels.py.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "assembly_analysis").is_dir()
)
sys.path.insert(0, str(Path(__file__).resolve().parent))

from quast_single_panels import render_shared_legend

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_quast_legend.png"


def main() -> int:
    render_shared_legend(OUTPUT_PNG)
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PNG.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
