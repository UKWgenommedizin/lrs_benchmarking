# Side analysis: does VACmap turn minimap2-unmapped reads into MAPQ 0 reads?

Exploratory hypothesis test. **Not part of the canonical benchmark.** Nothing here
feeds `tables/30x/final/` or `figures/30x/final/`, and no original CRAM/BAM is modified.

Started: 2026-09-30 · Status: **blocked on CRAM access** (preflight done, scripts ready)

## Hypothesis under test

> VACmap recovers difficult reads that minimap2 leaves unmapped; many originate from
> repetitive / low-mappability / structurally complex regions, and VACmap assigns them
> MAPQ 0 because their placement stays ambiguous.

The design must be able to reject this. The premise has two parts: (1) VACmap has fewer
unmapped reads, and (2) VACmap has more MAPQ 0 reads.

## Matched file pairs

See [`file_pairs.tsv`](file_pairs.tsv). Only these six comparisons are made. Samples and
technologies are never mixed.

| Sample | Tech | minimap2 (preset) | VACmap |
|---|---|---|---|
| HG002 | ONT | `HG002_ont_30x.hg38.mm2-ont.cram` (map-ont) | `HG002.ont.30x.hg38.vacmap-ont.cram` |
| HG003 | ONT | `HG003_ont_30x.hg38.mm2-ont.cram` (map-ont) | `HG003.ont.30x.hg38.vacmap-ont.cram` |
| HG004 | ONT | `HG004_ont_30x.hg38.mm2-ont.cram` (map-ont) | `HG004.ont.30x.hg38.vacmap-ont.cram` |
| HG002 | HiFi | `HG002_pb_30x.hg38.mm2-pb.cram` (map-hifi) | `HG002.pb.30x.hg38.vacmap-pb.cram` |
| HG003 | HiFi | `HG003_pb_30x.hg38.mm2-pb.cram` (map-hifi) | `HG003.pb.30x.hg38.vacmap-pb.cram` |
| HG004 | HiFi | `HG004_pb_30x.hg38.mm2-pb.cram` (map-hifi) | `HG004.pb.30x.hg38.vacmap-pb.cram` |

These pairs follow the presets used in the canonical benchmark table. **None of these CRAMs
exist on the local machine.** Paths come from the `samtools stats` headers:

- VACmap: `/data/genmedbfx/schilling_m/repos/lrs_benchmarking/cram/`
- minimap2: `/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/cram/`. The exact
  filenames are inferred and must be confirmed on the server. The PacBio stats headers
  record `HG00X_pb_30x.hg38.cram` for both minimap2 and pbmm2.

## Preflight result (from samtools stats, no CRAMs needed)

`scripts/00_stats_preflight.py` → `05_tables/00_stats_preflight_unmapped_accounting.tsv`

`reads unmapped` is **0 in every VACmap CRAM** because VACmap writes no record for reads it
cannot align. (`raw total sequences` counts primary records only, so this is not a
supplementary artefact.) Its unaligned reads are therefore `input reads − VACmap primary
records`:

| Sample | Tech | Input reads | mm2 unmapped | VACmap absent (unaligned) | ratio VAC/mm2 | mm2 MQ0 | VACmap MQ0 |
|---|---|---:|---:|---:|---:|---:|---:|
| HG002 | ONT | 5,203,769 | 85,034 (1.63%) | 143,911 (2.77%) | 1.69× | 22,248 | 86,635 |
| HG003 | ONT | 6,938,685 | 404,265 (5.83%) | 510,060 (7.35%) | 1.26× | 62,459 | 127,832 |
| HG004 | ONT | 7,046,468 | 304,942 (4.33%) | 415,653 (5.90%) | 1.36× | 58,199 | 111,031 |
| HG002 | HiFi | 6,736,453 | 16,859 (0.25%) | 35,393 (0.53%) | 2.10× | 10,300 | 86,715 |
| HG003 | HiFi | 6,621,583 | 14,247 (0.22%) | 30,700 (0.46%) | 2.15× | 9,259 | 81,763 |
| HG004 | HiFi | 6,491,546 | 13,829 (0.21%) | 25,967 (0.40%) | 1.88× | 9,167 | 65,151 |

**Observation:** once the absent reads are counted, VACmap leaves **more** reads unaligned
than minimap2 in all six datasets (1.3–2.2×). Premise (1) does not hold. It only appeared
to hold because VACmap reports 0 unmapped. Figure 11 already corrects for this. Premise (2)
holds: VACmap has 2–8× more MAPQ 0 reads.

**Consequence for the hypothesis:** VACmap cannot be lowering its unmapped rate, because it
isn't lower. The testable question becomes narrower: of the reads minimap2 leaves unmapped,
does VACmap place a subset at MAPQ 0, and where do those reads land? Stats cannot answer
this; read-level data is required. Upper bounds on the overlap (min of the two sets) are in
the TSV. In HiFi, VACmap MQ0 is 4–6× minimap2's unmapped count, so **most VACmap MQ0 reads
cannot come from minimap2-unmapped reads**, whatever the read-level overlap turns out to be.

*Caveat:* this assumes both aligners got the same FASTQ. minimap2, pbmm2 and VG Giraffe
agree exactly on input counts; VACmap was run on a different server. The read-level step
checks this directly (`reads_only_in_vacmap` must be 0).

## Validation issues to resolve on the server

1. **Reference build.** VACmap CRAMs were decoded against
   `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta`.
   minimap2 stats used `/mnt/storage/db/VarCAD_db/hg38/.../hg38.fasta.gz`. If the
   alignments used different references, MAPQ 0 is not comparable. The GIABv3 masking of
   false duplications lowers MQ0, and ALT contigs raise it. Check with:
   `samtools view -H X.cram | grep '^@SQ' | cut -f2,3,4 | md5sum`, compare `M5:` tags, and
   check `@PG` for the aligner command line.
2. **Same input reads.** Check that `reads_only_in_vacmap == 0` and that the minimap2 total
   equals the FASTQ read count.
3. **One primary per read.** `01_read_status.sh` logs `duplicate_qnames`, which must be 0.
4. **What "mapped" means for VACmap.** In the local test BAM, a VACmap primary record
   aligned 236 bp of a 28.6 kb read (99% clipped). `02_cross_tab.py` reports aligned
   fraction and clipping so that rescued reads with tiny aligned fragments are not counted
   as real rescues.

## Pipeline

| Step | Script | Status |
|---|---|---|
| 0 | `scripts/00_stats_preflight.py` | done (local) |
| 1 | `scripts/01_read_status.sh` - one streaming pass per CRAM, primary records only (`-F 0x900`, keeps unmapped) | written, smoke-tested on local 1k BAMs |
| 2 | `scripts/02_cross_tab.py` - cross-tab, forward and reciprocal outcomes, rescued-read MAPQ histogram, group BEDs (A–E), read characteristics, summary row | written, smoke-tested |
| 3 | repeat / difficult-region enrichment (bedtools + Fisher / OR with CI) | not written: needs annotations and bedtools |
| 4 | aggregation, figures, interpretation | after steps 1–3 |

Filtering choice: the proposed `-f 4 -F 2304` / `-F 2308` masks are correct. Step 1
instead uses a single `-F 0x900` pass that keeps primary mapped **and** primary unmapped
records. This yields the same read sets with one decode per CRAM instead of two.

### Commands (run on the server, from the repository root)

Estimated cost: one full decode of each of the 12 CRAMs (~30–60 min each at 8 threads)
plus a sort of about 7 M lines each (~1 GB RAM with `-S 3G`). Output is roughly
150–300 MB of gzip per CRAM.

```bash
D=alignment_analysis/side_analyses/vacmap_mapq0_hypothesis
REF=<reference used to write the CRAM>   # see validation issue 1

# 1a. sanity subset first (indexed region + unplaced unmapped)
bash $D/scripts/01_read_status.sh <mm2.cram> $REF $D/02_mapping_status/HG002.ONT.mm2.subset.tsv.gz 8 chr20:30000000-31000000 '*'
bash $D/scripts/01_read_status.sh <vac.cram> $REF $D/02_mapping_status/HG002.ONT.vac.subset.tsv.gz 8 chr20:30000000-31000000
python3 $D/scripts/02_cross_tab.py --sample HG002 --tech ONT_subset \
    --mm2 $D/02_mapping_status/HG002.ONT.mm2.subset.tsv.gz --vac $D/02_mapping_status/HG002.ONT.vac.subset.tsv.gz --outdir $D/07_logs/subset_check

# 1b. full run, per pair (6 pairs x 2 CRAMs)
bash $D/scripts/01_read_status.sh <mm2.cram> $REF $D/02_mapping_status/HG002.ONT.mm2.tsv.gz 8
bash $D/scripts/01_read_status.sh <vac.cram> $REF $D/02_mapping_status/HG002.ONT.vac.tsv.gz 8
python3 $D/scripts/02_cross_tab.py --sample HG002 --tech ONT \
    --mm2 $D/02_mapping_status/HG002.ONT.mm2.tsv.gz --vac $D/02_mapping_status/HG002.ONT.vac.tsv.gz --outdir $D
```

In the subset, reads whose other alignment falls outside the region show up as
`vac_absent`. It checks the logic and name compatibility, not the rates.

Only the per-read `*.tsv.gz` files (≈3 GB total) need to come back to this machine. The
CRAMs do not.

## Annotations needed (not available locally, nothing downloaded yet)

For GRCh38 with `chr` names, matching the GIABv3 analysis set:

| Priority | Annotation | Suggested source |
|---|---|---|
| 1 | RepeatMasker with class/family (LINE, SINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity) | UCSC hg38 `rmsk.txt.gz` |
| 2 | Segmental duplications | GIAB stratifications v3.x `GRCh38_segdups.bed.gz` (or UCSC `genomicSuperDups`) |
| 3 | Low mappability | GIAB `GRCh38_lowmappabilityall.bed.gz` |
| 4 | Centromeres / satellites | UCSC hg38 `centromeres` table (satellites also come from rmsk) |
| 5 | Tandem repeats / homopolymers | GIAB `GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz` |
| 6–7 | GIAB all-difficult union | GIAB `GRCh38_alldifficultregions.bed.gz` |
| 10 | CMRG | GIAB HG002 CMRG v1.00 GRCh38 benchmark BEDs |
| 10 | SV breakpoints | GIAB HG002 GRCh38 SV truth set (HG002 only, so no matching truth for HG003/HG004) |

Also required: `bedtools` (not installed in any local env).

## Interpretation rules

- Keep OBSERVATION, INTERPRETATION, HYPOTHESIS and EVIDENCE separate in every write-up.
- Report enrichment as proportions and odds ratios with CI. P-values alone are
  meaningless at n ≈ 10⁵.
- A MAPQ 0 read's VACmap coordinate is one arbitrary copy among equivalent loci. The
  repeat *class* is informative; the exact locus is not.
- Do not conclude "minimap2 discards repetitive reads". At most: "VACmap produced MAPQ 0
  alignments for a subset of minimap2-unmapped reads, enriched in …".
