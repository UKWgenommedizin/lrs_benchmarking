"""Read-level minimap2 vs VACmap comparison for one sample/technology pair.

Inputs are the per-read tables written by 01_read_status.sh. VACmap writes no
unmapped records, so a read present in the minimap2 table but absent from the
VACmap table is classified "vac_absent" (unaligned by VACmap). Reads present
only in VACmap mean the two runs did not use the same FASTQ -> reported and
the pair is flagged as failing validation.

Usage:
  02_cross_tab.py --sample HG002 --tech ONT --mm2 mm2.tsv.gz --vac vac.tsv.gz --outdir <analysis dir>

Writes (under <outdir>):
  02_mapping_status/{S}.{T}.crosstab.tsv          full minimap2 x VACmap status table
  02_mapping_status/{S}.{T}.mm2_unmapped_to_vac.tsv
  02_mapping_status/{S}.{T}.vac_unaligned_to_mm2.tsv   reciprocal
  02_mapping_status/{S}.{T}.rescued_vac_mapq_hist.tsv   VACmap MAPQ (all values) for mm2-unmapped reads
  01_read_lists/{S}.{T}.{group}.bed.gz            VACmap primary coordinates per analysis group
  04_statistics/{S}.{T}.read_characteristics.tsv
  05_tables/{S}.{T}.summary_row.tsv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

COLS = ["qname", "status", "mapq", "chrom", "start0", "end0",
        "read_len", "aligned_qlen", "clipped_bases", "nm", "n_supp"]
BINS = ["unmapped", "MQ0", "MQ1-9", "MQ10-19", "MQ20-29", "MQ30+"]
CONTROL_SAMPLE = 200_000  # confident-both control is subsampled for BED size


def load(path):
    d = pd.read_csv(path, sep="\t", names=COLS, na_values="NA",
                    dtype={"chrom": "category", "status": "category"})
    if d["qname"].duplicated().any():
        raise SystemExit(f"{path}: duplicate QNAMEs among primary records")
    return d.set_index("qname")


def mapq_bin(d, absent_label):
    b = pd.cut(d["mapq"], [-1, 0, 9, 19, 29, 1000], labels=BINS[1:]).astype(object)
    b[d["status"] == "unmapped"] = "unmapped"
    b[d["status"].isna()] = absent_label
    return b


def outcome_table(bins, order):
    n = bins.value_counts().reindex(order, fill_value=0)
    return pd.DataFrame({"outcome": order, "reads": n.values,
                         "pct": 100 * n.values / max(n.sum(), 1)})


def write_bed(d, path):
    bed = d.loc[d["status"] == "mapped", ["chrom", "start0", "end0", "mapq"]].copy()
    bed[["start0", "end0"]] = bed[["start0", "end0"]].astype(int)
    bed = bed.reset_index()[["chrom", "start0", "end0", "qname", "mapq"]]
    bed.sort_values(["chrom", "start0"]).to_csv(path, sep="\t", header=False, index=False)
    return len(bed)


def characteristics(d, group):
    m = d[d["status"] == "mapped"]
    aq = m["aligned_qlen"].replace(0, np.nan)
    return {
        "group": group, "reads": len(d),
        "median_read_len": d["read_len"].median(),
        "median_aligned_frac": (m["aligned_qlen"] / m["read_len"]).median(),
        "median_clip_frac": (m["clipped_bases"] / m["read_len"]).median(),
        "median_nm_per_kb_aligned": (1000 * m["nm"] / aq).median(),
        "pct_with_SA_tag": 100 * (m["n_supp"] > 0).mean() if len(m) else np.nan,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tech", required=True)
    ap.add_argument("--mm2", required=True)
    ap.add_argument("--vac", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    tag = f"{a.sample}.{a.tech}"
    for sub in ["01_read_lists", "02_mapping_status", "04_statistics", "05_tables"]:
        (out / sub).mkdir(parents=True, exist_ok=True)

    mm2, vac = load(a.mm2), load(a.vac)
    only_vac = vac.index.difference(mm2.index)
    if len(only_vac):
        print(f"VALIDATION FAIL: {len(only_vac)} VACmap reads absent from minimap2 "
              f"(different input FASTQ?) e.g. {list(only_vac[:3])}")

    j = mm2.join(vac, how="left", lsuffix="_mm2", rsuffix="_vac")
    split = lambda sfx: j[[c for c in j.columns if c.endswith(sfx)]].rename(columns=lambda c: c[: -len(sfx)])
    m, v = split("_mm2"), split("_vac")
    mb = mapq_bin(m, "absent")
    vb = mapq_bin(v, "unmapped")  # absent from VACmap output == unaligned by VACmap
    vb_detail = mapq_bin(v, "vac_absent")

    ct = pd.crosstab(mb, vb_detail).reindex(index=BINS, columns=["vac_absent"] + BINS, fill_value=0)
    ct.index.name, ct.columns.name = "minimap2", "VACmap"
    ct.to_csv(out / f"02_mapping_status/{tag}.crosstab.tsv", sep="\t")

    mm2_unm = mb == "unmapped"
    outcome_table(vb[mm2_unm], BINS).to_csv(out / f"02_mapping_status/{tag}.mm2_unmapped_to_vac.tsv", sep="\t", index=False)
    vac_unal = vb == "unmapped"
    outcome_table(mb[vac_unal], BINS).to_csv(out / f"02_mapping_status/{tag}.vac_unaligned_to_mm2.tsv", sep="\t", index=False)

    rescued = mm2_unm & (vb != "unmapped")
    v.loc[rescued, "mapq"].astype(int).value_counts().sort_index().rename_axis("vacmap_mapq").rename("reads") \
        .to_csv(out / f"02_mapping_status/{tag}.rescued_vac_mapq_hist.tsv", sep="\t")

    groups = {
        "A_mm2unmapped_vacMQ0": mm2_unm & (vb == "MQ0"),
        "C_mm2unmapped_vacMQ30plus": mm2_unm & (vb == "MQ30+"),
        "B_confident_both_MQ30plus": (mb == "MQ30+") & (vb == "MQ30+"),
        "D_mm2mapped_vacMQ0": (m["status"] == "mapped") & (vb == "MQ0"),
        "E_all_vac_mapped": vb != "unmapped",
    }
    rows = []
    for g, mask in groups.items():
        sub = v[mask]
        if g in ("B_confident_both_MQ30plus", "E_all_vac_mapped") and len(sub) > CONTROL_SAMPLE:
            sub = sub.sample(CONTROL_SAMPLE, random_state=1)
        write_bed(sub, out / f"01_read_lists/{tag}.{g}.bed.gz")
        rows.append(characteristics(sub, g))
    pd.DataFrame(rows).to_csv(out / f"04_statistics/{tag}.read_characteristics.tsv", sep="\t", index=False, float_format="%.4f")

    n_unm = int(mm2_unm.sum())
    resc = int(rescued.sum())
    vbu = vb[mm2_unm]
    summ = {
        "sample": a.sample, "technology": a.tech,
        "minimap2_total_reads": len(mm2),
        "minimap2_unmapped": n_unm,
        "vacmap_unaligned": int(vac_unal.sum()),
        "vacmap_mapq0": int((vb == "MQ0").sum()),
        "reads_only_in_vacmap": len(only_vac),
        "mm2_unmapped_vacmap_mapped": resc,
        "mm2_unmapped_vacmap_mapq0": int((vbu == "MQ0").sum()),
        "mm2_unmapped_vacmap_mapq10plus": int(vbu.isin(["MQ10-19", "MQ20-29", "MQ30+"]).sum()),
        "mm2_unmapped_vacmap_mapq30plus": int((vbu == "MQ30+").sum()),
        "pct_mm2_unmapped_rescued": 100 * resc / max(n_unm, 1),
        "pct_mm2_unmapped_to_vacmap_mapq0": 100 * (vbu == "MQ0").sum() / max(n_unm, 1),
        "pct_rescued_that_are_mapq0": 100 * (vbu == "MQ0").sum() / max(resc, 1),
        "pct_vacmap_mapq0_explained_by_mm2_unmapped": 100 * (vbu == "MQ0").sum() / max((vb == "MQ0").sum(), 1),
    }
    pd.DataFrame([summ]).to_csv(out / f"05_tables/{tag}.summary_row.tsv", sep="\t", index=False, float_format="%.4f")
    print(ct.to_string())
    print(pd.Series(summ).to_string())


if __name__ == "__main__":
    main()
