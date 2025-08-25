#!/usr/bin/env python3
# -*- coding: utf-8 -*-


##################
# IMPORT MODULES #
##################

import sys
import argparse
import logging
import signal
import re
import pysam
from pathlib import Path
from collections import defaultdict, Counter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Ignore SIGPIPE and handle it quietly
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

#############
# FUNCTIONS #
#############

def is_interactive():
    """Check if we're in an interactive session (e.g. Jupyter)."""
    import __main__ as main
    return not hasattr(main, '__file__')


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Convert GTF or BED12 to STAR SJ.out.tab-style junctions.\n"
            "Motif (col 5) is computed from a reference FASTA using pysam/htslib (.fai indexed). "
            "Outputs a failures TSV for any junctions whose bases could not be read.\n\n"
            "STAR column meanings used here:\n"
            " 1 chrom\n"
            " 2 intron start (1-based)\n"
            " 3 intron end (1-based)\n"
            " 4 strand: 0=undef, 1='+', 2='-'\n"
            " 5 motif on genomic + strand: 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT, 0 non-canonical\n"
            " 6 annotated (here: always 1)\n"
            " 7 uniquely mapping reads spanning junction (set via --uniq-reads)\n"
            " 8 multi-mapping reads spanning junction (here: 0)\n"
            " 9 maximum spliced overhang (here: 0)"
        )
    )
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--gtf',    help='Input GTF with exon features.')
    g.add_argument('--bed12',  help='Input BED12 with block structure (exons).')

    parser.add_argument('-o', '--out',        required=True,
                        help='Output SJ.out.tab path.')
    parser.add_argument('--ref-fasta',        required=True,
                        help='Reference FASTA (requires .fai index).')
    parser.add_argument('--fail',             default=None,
                        help='Failures TSV path (default: <out>.failures.tsv).')

    parser.add_argument('--min-intron',       type=int, default=20,
                        help='Drop introns shorter than this (default: 20).')
    parser.add_argument('--max-intron',       type=int, default=1_000_000,
                        help='Drop introns longer than this (default: 1,000,000).')
    parser.add_argument('--strand-agnostic',  action='store_true',
                        help='Force STAR strand code 0 for all junctions.')

    parser.add_argument(
        '--on-missing',
        choices=['skip', 'zero'],
        default='skip',
        help="When motif lookup fails (missing chrom, out-of-bounds, ambiguous bases): "
             "'skip' drop the junction (default); 'zero' keep it with motif=0 (STAR non-canonical)."
    )

    parser.add_argument(
        '--uniq-reads', type=int, default=100,
        help="Value to write into SJ.out.tab column 7 (unique spanning reads). "
             "Some downstream tools ignore junctions with 0. Default: 100."
    )

    # Show help if no args when running from CLI; empty args in notebooks.
    if is_interactive():
        args = parser.parse_args('')
    else:
        args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args


def parse_gtf(path):
    """
    Yield (chrom, strand, [(exon_start_1b, exon_end_1b), ...]) per transcript.
    GTF exon coordinates are 1-based inclusive.
    """
    tx_exons = defaultdict(list)
    tx_strand = {}
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != 'exon':
                continue
            m = re.search(r'transcript_id\s+"([^"]+)"', attrs) or re.search(r"transcript_id\s+'([^']+)'", attrs)
            if not m:
                m2 = re.search(r'gene_id\s+"([^"]+)"', attrs) or re.search(r"gene_id\s+'([^']+)'", attrs)
                m3 = re.search(r'exon_number\s+"?(\d+)"?', attrs)
                if m2 and m3:
                    txid = f"{m2.group(1)}:exonset"
                else:
                    continue
            else:
                txid = m.group(1)
            s = int(start)
            e = int(end)
            tx_exons[(chrom, txid)].append((s, e))
            if (chrom, txid) not in tx_strand:
                tx_strand[(chrom, txid)] = strand
    for (chrom, txid), exons in tx_exons.items():
        yield chrom, tx_strand[(chrom, txid)], sorted(exons, key=lambda x: x[0])


def parse_bed12(path):
    """
    Yield (chrom, strand, [(exon_start_1b, exon_end_1b), ...]) per BED12 line.
    BED is 0-based, half-open; convert to 1-based inclusive for exons.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            chrom = parts[0]
            chrom_start = int(parts[1])
            strand = parts[5] if parts[5] in ('+','-') else '.'
            block_count = int(parts[9])
            block_sizes = [int(x) for x in parts[10].rstrip(',').split(',')]
            block_starts = [int(x) for x in parts[11].rstrip(',').split(',')]
            if len(block_sizes) != block_count or len(block_starts) != block_count:
                continue
            exons = []
            for bs, bstart in zip(block_sizes, block_starts):
                exon_start0 = chrom_start + bstart
                exon_end0   = exon_start0 + bs
                exon_start1 = exon_start0 + 1
                exon_end1   = exon_end0
                exons.append((exon_start1, exon_end1))
            yield chrom, strand, sorted(exons, key=lambda x: x[0])


def introns_from_exons(exons_1b):
    """
    Given sorted exons as (start,end) 1-based inclusive, return introns as (start,end) 1-based inclusive.
    Intron start = prev_exon_end + 1; intron end = next_exon_start - 1
    """
    introns = []
    for i in range(len(exons_1b) - 1):
        prev_end   = exons_1b[i][1]
        next_start = exons_1b[i+1][0]
        intr_start = prev_end + 1
        intr_end   = next_start - 1
        if intr_end >= intr_start:
            introns.append((intr_start, intr_end))
    return introns


def strand_to_code(s):
    if s == '+': return 1
    if s == '-': return 2
    return 0


def motif_code_from_dinucs(donor, acceptor):
    """
    STAR SJ.out.tab column 5 (on genomic + strand):
      1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT, 0 non-canonical
    """
    pair = (donor, acceptor)
    if pair == ('GT','AG'): return 1
    if pair == ('CT','AC'): return 2
    if pair == ('GC','AG'): return 3
    if pair == ('CT','GC'): return 4
    if pair == ('AT','AC'): return 5
    if pair == ('GT','AT'): return 6
    return 0  # non-canonical


def fetch_dinucs(ref, chrom, intr_start_1b, intr_end_1b):
    """
    Fetch donor and acceptor dinucleotides from indexed FASTA using pysam.
    Coordinates are 1-based inclusive. pysam.fetch uses 0-based, end-exclusive.
      donor    = [start, start+1] -> fetch(start-1, start+1)
      acceptor = [end-1, end]     -> fetch(end-2, end)
    Returns (donor, acceptor, reason_if_failed_or_None)
    """
    try:
        donor    = ref.fetch(chrom, intr_start_1b - 1, intr_start_1b + 1).upper()
        acceptor = ref.fetch(chrom, intr_end_1b - 2,   intr_end_1b      ).upper()
    except (KeyError, ValueError) as e:
        return (None, None, f"lookup_failed:{e}")
    if len(donor) != 2 or len(acceptor) != 2:
        return (None, None, f"length_issue:{chrom}:{intr_start_1b}-{intr_end_1b}")
    if any(b not in 'ACGT' for b in donor + acceptor):
        return (None, None, f"ambiguous_bases:{donor}/{acceptor}")
    return donor, acceptor, None


def write_sj(path_out, rows_iterable):
    with open(path_out, 'w') as out:
        for row in rows_iterable:
            out.write('\t'.join(map(str, row)) + '\n')


########
# MAIN #
########

def main():
    args = parse_args()

    # Input existence checks
    in_path = Path(args.gtf) if args.gtf else Path(args.bed12)
    if not in_path.is_file():
        logging.error(f"File not found: {in_path}")
        sys.exit(1)

    ref_path = Path(args.ref_fasta)
    if not ref_path.is_file():
        logging.error(f"Reference FASTA not found: {ref_path}")
        sys.exit(1)
    fai_path = Path(str(ref_path) + '.fai')
    if not fai_path.is_file():
        logging.info("FASTA index (.fai) not found; pysam/htslib will attempt to create it if permissions allow.")

    # Open FASTA
    logging.info(f"Opening reference FASTA: {ref_path}")
    try:
        ref = pysam.FastaFile(str(ref_path))
    except Exception as e:
        logging.error(f"Failed to open FASTA via pysam: {e}")
        sys.exit(1)

    # Iterator over transcripts
    iterator = parse_gtf(str(in_path)) if args.gtf else parse_bed12(str(in_path))

    # Failures file
    fail_path = args.fail if args.fail else (args.out + '.failures.tsv')
    try:
        fail_fh = open(fail_path, 'w')
    except Exception as e:
        logging.error(f"Could not open failures file for writing: {e}")
        sys.exit(1)
    fail_fh.write('chrom\tstart\tend\treason\n')

    sj_set = set()
    n_total_introns = 0
    n_kept = 0
    n_failed = 0
    n_skipped_by_len = 0
    n_skipped_missing = 0
    motif_counts = Counter()

    logging.info("Scanning transcripts and computing junction motifs...")
    for chrom, strand, exons in iterator:
        introns = introns_from_exons(exons)
        s_code = 0 if args.strand_agnostic else strand_to_code(strand)
        for intr_start, intr_end in introns:
            n_total_introns += 1
            length = intr_end - intr_start + 1
            if length < args.min_intron or length > args.max_intron:
                n_skipped_by_len += 1
                continue

            donor, acceptor, reason = fetch_dinucs(ref, chrom, intr_start, intr_end)
            if reason is None:
                motif_val = motif_code_from_dinucs(donor, acceptor)
            else:
                n_failed += 1
                fail_fh.write(f"{chrom}\t{intr_start}\t{intr_end}\t{reason}\n")
                if args.on_missing == 'skip':
                    n_skipped_missing += 1
                    continue
                motif_val = 0  # user explicitly chose to emit non-canonical (0) on missing

            # STAR SJ.out.tab columns:
            # 1:chrom 2:start 3:end 4:strand(0/1/2) 5:motif 6:annotated 7:uniq 8:multi 9:max_overhang
            row = (
                chrom, intr_start, intr_end,
                s_code, motif_val,
                1,                # annotated (present in annotation)
                args.uniq_reads,  # unique spanning reads
                0,                # multi-mapping spanning reads
                0                 # max overhang
            )
            sj_set.add(row)
            n_kept += 1
            motif_counts[motif_val] += 1

    # Write output
    logging.info(f"Writing SJ.out.tab to: {args.out}")
    write_sj(args.out, sorted(sj_set, key=lambda r: (r[0], r[1], r[2], r[3])))

    fail_fh.close()
    logging.info(f"Failures logged to: {fail_path}")

    # Motif summary
    canonical_total = sum(motif_counts[m] for m in (1,2,3,4,5,6))
    noncanonical_total = motif_counts[0]
    logging.info("Motif distribution (col 5): " +
                 ", ".join(f"{m}={motif_counts[m]}" for m in range(0,7)))
    logging.info(f"Canonical total (1–6): {canonical_total}; Non-canonical (0): {noncanonical_total}")

    logging.info(
        f"Done. introns_total={n_total_introns} kept={n_kept} "
        f"skipped_by_length={n_skipped_by_len} lookup_failures_logged={n_failed} "
        f"skipped_by_missing={n_skipped_missing}"
    )


if __name__ == '__main__':
    main()
