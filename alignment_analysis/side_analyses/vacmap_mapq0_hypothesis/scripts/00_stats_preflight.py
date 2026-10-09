"""Preflight read accounting for the VACmap-MAPQ0 hypothesis, from samtools stats only.

Uses the canonical benchmark table (samtools stats SN values, which count primary
reads only: "raw total sequences ... excluding supplementary and secondary reads").

VACmap writes no unmapped records, so its unaligned reads are recovered as
    shared raw input (minimap2 raw_total_sequences) - VACmap raw_total_sequences.
This assumes both aligners received the same FASTQ; the read-level step
(01_read_status.sh + 02_cross_tab.py) verifies that assumption by QNAME.

Output: 05_tables/00_stats_preflight_unmapped_accounting.tsv
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[2]
TABLE = REPO / "alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv"
OUT = HERE / "05_tables/00_stats_preflight_unmapped_accounting.tsv"

d = pd.read_csv(TABLE, sep="\t")
d["aligner"] = d["aligner"].replace({"VACMap": "VACmap"})
cols = ["sample", "read_technology", "raw_total_sequences", "reads_unmapped", "reads_mq0"]
mm2 = d.loc[d["aligner"] == "minimap2", cols].set_index(["sample", "read_technology"])
vac = d.loc[d["aligner"] == "VACmap", cols].set_index(["sample", "read_technology"])

out = pd.DataFrame(index=mm2.index)
out["input_reads"] = mm2["raw_total_sequences"]
out["mm2_unmapped"] = mm2["reads_unmapped"]
out["mm2_mq0"] = mm2["reads_mq0"]
out["vac_primary_records"] = vac["raw_total_sequences"]
out["vac_unmapped_reported"] = vac["reads_unmapped"]
out["vac_absent_from_output"] = out["input_reads"] - out["vac_primary_records"]
out["vac_mq0"] = vac["reads_mq0"]
out["pct_mm2_unmapped"] = 100 * out["mm2_unmapped"] / out["input_reads"]
out["pct_vac_unaligned"] = 100 * out["vac_absent_from_output"] / out["input_reads"]
out["vac_unaligned_over_mm2_unmapped"] = out["vac_absent_from_output"] / out["mm2_unmapped"]
out["pct_mm2_mq0_of_input"] = 100 * out["mm2_mq0"] / out["input_reads"]
out["pct_vac_mq0_of_input"] = 100 * out["vac_mq0"] / out["input_reads"]
# Upper bound on "mm2-unmapped -> VACmap MQ0": cannot exceed either set.
out["max_possible_overlap"] = out[["mm2_unmapped", "vac_mq0"]].min(axis=1)

out = out.reset_index().sort_values(["read_technology", "sample"])
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
print(out.to_string(index=False))
print(f"\nwrote {OUT.relative_to(REPO)}")
