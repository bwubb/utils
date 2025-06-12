#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

def parse_dict(dict_path):
    contig_lengths = {}
    with open(dict_path) as f:
        for line in f:
            if line.startswith("@SQ"):
                fields = dict(field.split(":", 1) for field in line.strip().split("\t")[1:] if ":" in field)
                if "SN" in fields and "LN" in fields:
                    contig_lengths[fields["SN"]] = int(fields["LN"])
    return contig_lengths

def detect_prefix_style(names):
    """Returns 'chr' if majority names start with chr, else ''."""
    chr_prefixed = sum(name.startswith("chr") for name in names)
    return "chr" if chr_prefixed >= len(names) / 2 else ""

def normalize_contig(name, desired_prefix):
    base = name[3:] if name.startswith("chr") else name
    return f"chr{base}" if desired_prefix == "chr" else base

def find_dict_for_fasta(fasta_path):
    dict_path = Path(str(fasta_path) + ".dict")
    if not dict_path.exists():
        raise FileNotFoundError(f"No .dict found for FASTA at {dict_path}")
    return dict_path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--intervals", required=True, help="Interval list file")
    parser.add_argument("-r", "--reference", required=True, help="FASTA file or .dict file")
    parser.add_argument("-o", "--output", required=True, help="Output file for cleaned intervals")
    args = parser.parse_args()

    dict_file = Path(args.reference)
    if dict_file.suffix != ".dict":
        dict_file = find_dict_for_fasta(dict_file)

    contig_lengths = parse_dict(dict_file)

    with open(args.intervals) as fin:
        interval_lines = [line.strip() for line in fin if re.match(r"^\S+:\d+-\d+$", line)]

    interval_chroms = {line.split(":")[0] for line in interval_lines}
    contig_chroms = set(contig_lengths.keys())

    interval_style = detect_prefix_style(interval_chroms)
    dict_style = detect_prefix_style(contig_chroms)

    # Normalize contig lengths to interval's naming style
    norm_contig_lengths = {
        normalize_contig(name, interval_style): length
        for name, length in contig_lengths.items()
    }

    kept, removed = 0, 0
    with open(args.output, "w") as fout:
        for line in interval_lines:
            chrom, pos = line.split(":")
            start, end = map(int, pos.split("-"))
            if chrom in norm_contig_lengths and end <= norm_contig_lengths[chrom]:
                fout.write(line + "\n")
                kept += 1
            else:
                removed += 1

    print(f"[INFO] Wrote {kept} valid intervals to: {args.output}")
    print(f"[INFO] Skipped {removed} invalid interval(s) that exceeded contig bounds.")

if __name__ == "__main__":
    main()
