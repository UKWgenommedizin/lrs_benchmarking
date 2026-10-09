#!/usr/bin/env bash
# One streaming pass over a CRAM/BAM -> one line per read (primary record only).
#
# Usage: 01_read_status.sh <in.cram> <reference.fasta> <out.tsv.gz> [threads] [regions...]
#   regions (optional, needs .crai/.bai): e.g. chr20:30000000-31000000 '*'
#   for the subset sanity check; '*' = unplaced unmapped reads. Omit for the whole file.
#
# Filter: -F 0x900 drops secondary (0x100) and supplementary (0x800) records and
# KEEPS primary unmapped records (0x4), so a single pass yields both mapped and
# unmapped reads. Every read has exactly one primary record in a valid SAM; the
# duplicate-QNAME check below verifies that for each aligner (VACmap included).
#
# Output columns (sorted by qname, LC_ALL=C):
#   qname status mapq chrom start0 end0 read_len aligned_qlen clipped_bases nm n_supp_tag
#   status: unmapped | mapped.  read_len counts H-clipped bases too.
#   start0/end0 are 0-based half-open reference coordinates (BED-ready).
set -euo pipefail

in=$1; ref=$2; out=$3; threads=${4:-8}; shift $(( $# < 4 ? $# : 4 )); regions=("$@")
tmp=$(dirname "$out")/tmp_sort; mkdir -p "$tmp"
log=${out%.tsv.gz}.log

samtools view -@ "$threads" -T "$ref" -F 0x900 "$in" "${regions[@]}" \
| awk -F'\t' -v OFS='\t' '
  {
    unm = int($2 / 4) % 2
    rl = 0; aq = 0; clip = 0; rlen = 0; cig = $6
    while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
      n = substr(cig, 1, RLENGTH - 1) + 0; op = substr(cig, RLENGTH, 1)
      if (op ~ /[MI=X]/) aq += n
      if (op ~ /[MI=XS]/) rl += n
      if (op == "H") rl += n
      if (op ~ /[SH]/) clip += n
      if (op ~ /[MDN=X]/) rlen += n
      cig = substr(cig, RLENGTH + 1)
    }
    if ($6 == "*") rl = length($10)
    nm = "NA"; sa = 0
    for (i = 12; i <= NF; i++) {
      if (substr($i, 1, 5) == "NM:i:") nm = substr($i, 6)
      if (substr($i, 1, 5) == "SA:Z:") sa = gsub(/;/, ";", $i)
    }
    if (unm) print $1, "unmapped", "NA", "NA", "NA", "NA", rl, 0, 0, "NA", 0
    else     print $1, "mapped", $5, $3, $4 - 1, $4 - 1 + rlen, rl, aq, clip, nm, sa
  }' \
| LC_ALL=C sort -t$'\t' -k1,1 -S 3G -T "$tmp" --parallel="$threads" \
| gzip -1 > "$out"

# Validation: every QNAME must appear once.
{
  echo "input	$in"
  echo "primary_records	$(zcat "$out" | wc -l)"
  echo "unique_qnames	$(zcat "$out" | cut -f1 | LC_ALL=C uniq | wc -l)"
  echo "duplicate_qnames	$(zcat "$out" | cut -f1 | LC_ALL=C uniq -d | wc -l)"
  echo "unmapped	$(zcat "$out" | awk '$2=="unmapped"' | wc -l)"
  echo "mapq0	$(zcat "$out" | awk '$2=="mapped" && $3==0' | wc -l)"
} | tee "$log"
