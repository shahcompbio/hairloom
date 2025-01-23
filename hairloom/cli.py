import pysam
import pandas as pd
import argparse

def cli():
    parser = argparse.ArgumentParser(
        description="Hairloom: Analyze split reads of long-read sequencing data."
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # Define the `extract` subcommand
    extract_parser = subparsers.add_parser("extract", help="Extract split-read alignments")
    extract_parser.add_argument("bam", type=str, help="Input BAM file")
    extract_parser.add_argument("chrom", type=str, help="Chromosome name (e.g., chr1)")
    extract_parser.add_argument("start", type=int, help="Start position (1-based)")
    extract_parser.add_argument("end", type=int, help="End position (1-based)")

    # `breakpoints` subcommand
    breakpoints_parser = subparsers.add_parser("breakpoints", help="Generate breakpoint table")
    breakpoints_parser.add_argument("bam", type=str, help="Input BAM file")
    breakpoints_parser.add_argument("chrom", type=str, help="Chromosome name (e.g., chr1)")
    breakpoints_parser.add_argument("start", type=int, help="Start position (1-based)")
    breakpoints_parser.add_argument("end", type=int, help="End position (1-based)")

    # `segments` command
    segments_parser = subparsers.add_parser("segments", help="Generate segment table")
    segments_parser.add_argument("bam", type=str, help="Input BAM file")
    segments_parser.add_argument("chrom", type=str, help="Chromosome name (e.g., chr1)")
    segments_parser.add_argument("start", type=int, help="Start position (1-based)")
    segments_parser.add_argument("end", type=int, help="End position (1-based)")

    # `svs` command
    svs_parser = subparsers.add_parser("svs", help="Generate SV tables for paired breakpoints")
    svs_parser.add_argument("bam", type=str, help="Input BAM file")
    svs_parser.add_argument("chrom", type=str, help="Chromosome name (e.g., chr1)")
    svs_parser.add_argument("start", type=int, help="Start position (1-based)")
    svs_parser.add_argument("end", type=int, help="End position (1-based)")

    args = parser.parse_args()

    if args.command == "extract":
        handle_extract(args)
    elif args.command == "breakpoints":
        handle_breakpoints(args)
    elif args.command == "segments":
        handle_segments(args)
    elif args.command == "svs":
        handle_svs(args)


def handle_extract(args):
    """Handle the 'extract' command."""
    try:
        bam_file = pysam.AlignmentFile(args.bam, "rb")
    except FileNotFoundError:
        print(f"Error: BAM file '{args.bam}' not found.")
        return

    from hairloom import extract_read_data
    read_table = extract_read_data(bam_file, args.chrom, args.start, args.end)

    # Output the read table as a TSV
    print(read_table.to_csv(sep="\t", index=False))


def handle_breakpoints(args):
    """Handle the 'breakpoints' command."""
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    from hairloom import extract_read_data, make_bundle
    read_table = extract_read_data(bam_file, args.chrom, args.start, args.end)
    bundle = make_bundle(read_table)

    # Generate breakpoint table
    from hairloom import make_brk_table
    table = make_brk_table(bundle)

    # Output the breakpoint table as TSV
    print(table.to_csv(sep="\t", index=False))


def handle_segments(args):
    """Handle the 'extract' command."""
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    from hairloom import extract_read_data, make_bundle
    read_table = extract_read_data(bam_file, args.chrom, args.start, args.end)
    bundle = make_bundle(read_table)

    # Generate breakpoint table
    from hairloom import make_seg_table
    table = make_seg_table(bundle)

    # Output the breakpoint table as TSV
    print(table.to_csv(sep="\t", index=False))


def handle_svs(args):
    """Handle the 'extract' command."""
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    from hairloom import extract_read_data, make_bundle
    read_table = extract_read_data(bam_file, args.chrom, args.start, args.end)
    bundle = make_bundle(read_table)

    # Generate breakpoint table
    from hairloom import make_tra_table
    table = make_tra_table(bundle)

    # Output the breakpoint table as TSV
    print(table.to_csv(sep="\t", index=False))