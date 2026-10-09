#!/usr/bin/env python3
"""Extract every table from the F2 report and check whether its numbers exist in repository data.

For each captioned ``table`` environment in the LaTeX source, the script:
  1. parses the ``tabular`` body into rows,
  2. searches tabular data files in the repository (tsv/csv/txt, excluding this script's
     own output directory) for each row,
     requiring all numeric cells of the row to appear as exact tokens on one line,
  3. writes the extracted table and a provenance verdict to ``final_report_files/tables/``.

A row is VERIFIED only if a data file (not a LaTeX file or a LaTeX build log) contains it.
Outputs:
  tables/report_table<N>_extracted.tsv   the table exactly as printed in the report
  tables/report_tables_provenance.tsv    one line per row with verdict and matching files
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_TEX = REPO / "final_report_files" / "final_f2_report.tex"
DEFAULT_OUT = REPO / "final_report_files" / "tables"

DATA_SUFFIXES = {".tsv", ".csv", ".txt"}
SKIP_DIRS = {".git", ".snakemake", "__pycache__", "node_modules", "conda", "tutorials", "tables_out"}
# Files that only echo LaTeX source (the report itself, the command cheat-sheet, pdflatex logs).
NON_DATA_PARTS = ("build_logs",)
MAX_BYTES = 200 * 1024 * 1024
TOKEN_SPLIT = re.compile(r"[\s,;|]+")


def strip_latex(cell):
    cell = re.sub(r"\\(textbf|textit|emph|texttt|mathrm)\{([^}]*)\}", r"\2", cell)
    cell = re.sub(r"\\(toprule|midrule|bottomrule|hline)", "", cell)
    return cell.replace("$", "").replace("{", "").replace("}", "").strip()


def extract_tables(tex_text):
    """Yield (index, caption, label, header, rows) for each table environment."""
    body = re.sub(r"(?<!\\)%.*", "", tex_text)
    for n, env in enumerate(re.findall(r"\\begin\{table\}(.*?)\\end\{table\}", body, re.S), start=1):
        tab = re.search(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}", env, re.S)
        if not tab:
            continue
        caption = re.search(r"\\caption\{(.*?)\}\s*(\\label|$)", env, re.S)
        label = re.search(r"\\label\{(.*?)\}", env)
        lines = [ln for ln in re.split(r"\\\\", tab.group(1))]
        rows = []
        for ln in lines:
            cells = [strip_latex(c) for c in ln.split("&")]
            if any(cells):
                rows.append(cells)
        yield (n, caption.group(1).strip() if caption else "",
               label.group(1) if label else "", rows[0], rows[1:])


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def candidate_files(root, exclude):
    for dirpath, dirnames, filenames in os.walk(root):
        if Path(dirpath).resolve() == exclude:
            dirnames[:] = []
            continue
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.endswith("_")]
        for name in filenames:
            p = Path(dirpath, name)
            if p.suffix.lower() in DATA_SUFFIXES:
                try:
                    if p.stat().st_size <= MAX_BYTES:
                        yield p
                except OSError:
                    continue


def search_rows(rows, root, exclude):
    """Return {row_index: [(relative_path, line_no), ...]} for rows whose numeric cells co-occur on one line."""
    targets = []
    for i, row in enumerate(rows):
        nums = {c for c in row if is_number(c)}
        if nums:
            # longest value is the cheapest substring pre-filter before tokenising the line
            targets.append((i, nums, max(nums, key=len)))
    hits = {i: [] for i, _, _ in targets}
    for path in candidate_files(root, exclude):
        try:
            with open(path, "r", errors="ignore") as fh:
                for line_no, line in enumerate(fh, start=1):
                    tokens = None
                    for i, nums, probe in targets:
                        if probe not in line:
                            continue
                        tokens = tokens or set(TOKEN_SPLIT.split(line))
                        if nums <= tokens:
                            hits[i].append((str(path.relative_to(root)), line_no))
        except OSError:
            continue
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--root", type=Path, default=REPO)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    tables = list(extract_tables(args.tex.read_text()))
    if not tables:
        sys.exit(f"No table environments found in {args.tex}")

    prov_path = args.out / "report_tables_provenance.tsv"
    all_verified = True
    with open(prov_path, "w", newline="") as prov_fh:
        prov = csv.writer(prov_fh, delimiter="\t", lineterminator="\n")
        prov.writerow(["table", "label", "caption", "row", "values", "verdict", "data_file_matches", "non_data_matches"])
        for n, caption, label, header, rows in tables:
            with open(args.out / f"report_table{n}_extracted.tsv", "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(header)
                w.writerows(rows)

            hits = search_rows(rows, args.root, args.out.resolve())
            print(f"Table {n} ({label}): {caption}")
            for i, row in enumerate(rows):
                found = hits.get(i, [])
                data = [f"{p}:{ln}" for p, ln in found
                        if not p.endswith(".tex") and not any(part in p for part in NON_DATA_PARTS)]
                other = [f"{p}:{ln}" for p, ln in found if f"{p}:{ln}" not in data]
                verdict = "VERIFIED" if data else "NOT_FOUND_IN_DATA"
                all_verified &= bool(data)
                prov.writerow([n, label, caption, i + 1, " | ".join(row), verdict,
                               ";".join(data) or "-", ";".join(other) or "-"])
                print(f"  row {i + 1}: {' | '.join(row):40s} {verdict}"
                      + (f"  ({len(data)} data hits)" if data else ""))

    print(f"\nWrote {prov_path.relative_to(args.root)}")
    return 0 if all_verified else 1


if __name__ == "__main__":
    sys.exit(main())
