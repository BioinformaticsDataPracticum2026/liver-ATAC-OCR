#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import gzip

# =========================
# 1. config
# =========================
PROJECT_DIR = Path("/Users/yefenglin/Desktop/data_class/project")

HUMAN_INPUT = PROJECT_DIR / "human_liver.idr.optimal_peak.narrowPeak.gz"
MOUSE_INPUT = PROJECT_DIR / "mouse_liver.idr.optimal_peak.narrowPeak.gz"

OUTPUT_DIR = PROJECT_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

HUMAN_CLEAN_BED = OUTPUT_DIR / "human_liver.clean.bed"
MOUSE_CLEAN_BED = OUTPUT_DIR / "mouse_liver.clean.bed"

HUMAN_SUMMIT_BED = OUTPUT_DIR / "human_liver.summits.bed"
MOUSE_SUMMIT_BED = OUTPUT_DIR / "mouse_liver.summits.bed"


# =========================
# 2. utility functions
# =========================
def open_maybe_gzip(file_path: Path):
    """
    Automatically open either a plain text file or a gzip-compressed file.
    Even though your filename ends with .gz.gz, it will still be read as gzip.
    """
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def chr_sort_key(chrom: str):
    """
    Sort chromosomes in a more natural order:
    chr1, chr2, ..., chr22, chrX, chrY, chrM, others
    """
    if chrom.startswith("chr"):
        name = chrom[3:]
    else:
        name = chrom

    if name.isdigit():
        return (0, int(name))

    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name in special:
        return (1, special[name])

    return (2, name)


def sort_bed_records(records):
    """
    Sort BED records by chromosome, start, and end position.
    """
    return sorted(records, key=lambda x: (chr_sort_key(x[0]), x[1], x[2]))


def parse_narrowpeak(input_file: Path):
    """
    Read a narrowPeak file and return two lists:
    1. clean peaks:   (chr, start, end)
    2. summits:       (chr, summit, summit+1)
    """
    peaks = []
    summits = []

    with open_maybe_gzip(input_file) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()

            # Skip empty lines, comments, and header lines
            if not line:
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            fields = line.split("\t")
            if len(fields) < 10:
                print(f"[WARNING] {input_file.name} line {line_num} is not a standard narrowPeak record (fewer than 10 columns), skipped.")
                continue

            chrom = fields[0]

            # Remove chrM peaks
            if chrom == "chrM":
                continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                print(f"[WARNING] {input_file.name} line {line_num} has non-integer start/end coordinates, skipped.")
                continue

            peaks.append((chrom, start, end))

            # Column 10 in narrowPeak: summit offset relative to peak start
            try:
                summit_offset = int(fields[9])
            except ValueError:
                summit_offset = -1

            if summit_offset != -1:
                summit_pos = start + summit_offset
                if start <= summit_pos < end:
                    summits.append((chrom, summit_pos, summit_pos + 1))

    return peaks, summits


def write_bed(records, output_file: Path):
    """
    Write BED3 records to a file.
    """
    with open(output_file, "w") as f:
        for chrom, start, end in records:
            f.write(f"{chrom}\t{start}\t{end}\n")


def process_one(input_file: Path, clean_out: Path, summit_out: Path, label: str):
    """
    Process one species.
    """
    print(f"\n[INFO] Processing {label}")
    print(f"[INFO] Input file: {input_file}")

    if not input_file.exists():
        raise FileNotFoundError(f"{label} input file not found: {input_file}")

    peaks, summits = parse_narrowpeak(input_file)

    peaks = sort_bed_records(peaks)
    summits = sort_bed_records(summits)

    write_bed(peaks, clean_out)
    write_bed(summits, summit_out)

    print(f"[DONE] {label} finished")
    print(f"       clean peaks : {clean_out} ({len(peaks)} lines)")
    print(f"       summits     : {summit_out} ({len(summits)} lines)")


# =========================
# 3. main program
# =========================
def main():
    print("==== ATAC peak preprocessing start ====")
    print(f"Project dir: {PROJECT_DIR}")

    process_one(HUMAN_INPUT, HUMAN_CLEAN_BED, HUMAN_SUMMIT_BED, "human_liver")
    process_one(MOUSE_INPUT, MOUSE_CLEAN_BED, MOUSE_SUMMIT_BED, "mouse_liver")

    print("\n==== ALL DONE ====")
    print("Output files:")
    print(HUMAN_CLEAN_BED)
    print(HUMAN_SUMMIT_BED)
    print(MOUSE_CLEAN_BED)
    print(MOUSE_SUMMIT_BED)


if __name__ == "__main__":
    main()