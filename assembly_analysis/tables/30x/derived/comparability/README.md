# Assembler comparability audit (30x)

Audit of whether the QUAST values behind `assembler_benchmark_bars.png` are comparable across Flye, GoldRush and Verkko. All values come from the local copies of `assembly_quality/quast/{assembler}/{dataset}/report.tsv`. The assembly FASTAs, `quast.log`, Flye `assembly_info.txt` and assembler logs are only on the server, so FASTA-level values are NA until `fasta_comparability_stats.py` is run there.

| File | Content |
|---|---|
| `assembler_comparability_audit.tsv` | One row per sample x assembler x technology: FASTA path, representation, original QUAST values, comparability flags, evidence |
| `assembler_standardized_metrics.tsv` | Values eligible for the main figure; NA where comparability is not established |
| `quast_contig_redundancy_30x.tsv` | Reference depth and contained (haplotig-like) contigs per assembly, from QUAST `contigs_reports/all_alignments_assembly.tsv` (`scripts/metrics/quast_contig_redundancy.py`) |
| `fasta_comparability_stats.tsv` | *(produced on the server)* scaffold vs gap-split contig statistics, and Verkko header-prefix / haplotype statistics |
| `flye_assembly_info_summary.tsv` | *(produced on the server)* Flye `alt_group` and coverage summary for the redundancy check |

## Findings

- **Flye ONT:** haploid-like (about 2.95 Gb, duplication 1.01, no N gaps within QUAST rounding). Its contig N50 and count are comparable.
- **Flye HiFi:** 3.95–3.97 Gb and duplication 1.33–1.34. The workflow (`--pacbio-hifi`, without `--keep-haplotypes`) should give a collapsed assembly. *Possible retained redundant/haplotig sequence — not yet confirmed.* Contiguity and the absolute misassembly counts stay NA until `assembly_info.txt` has been checked.
- **GoldRush:** the final FASTA is the ntLink 5-round scaffold output, with 157–233 N per 100 kbp. QUAST ran without `--split-scaffolds`, so its N50 and count are **scaffold** values. Contig values need the gap split. The reference metrics are comparable: QUAST reports scaffold-gap misassemblies separately.
- **Verkko:** run without Hi-C, trio or Pore-C input, so the workflow keeps only the combined `assembly.fasta` (about 6.1 Gb, duplication 2.04–2.09). QUAST compared this **combined diploid** assembly with haploid GRCh38. No per-haplotype FASTAs exist, and pooled N50 is not comparable. Misassembly counts cover about twice the sequence.

## To complete the audit (on the server, from the repository root)

```bash
python assembly_analysis/scripts/metrics/fasta_comparability_stats.py
```

This writes the two server-produced tables listed above. It also writes the GoldRush gap-split FASTAs to `assembly_quality/derived/comparability/goldrush_gapsplit/{dataset}/assembly.gapsplit.fasta`. The original assemblies are only read.
