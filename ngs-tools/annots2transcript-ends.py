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
    desc = (
        "Extract transcription start sites (TSS) or transcription termination sites (TTS) "
        "from an annotation file (GTF or BED12) and write them in BED6 format.\n\n"
        "BED6 columns:\n"
        " 1 chrom\n 2 start (0-based)\n 3 end (0-based, half-open)\n 4 name\n 5 score (0-1000)\n 6 strand (+/-)\n\n"
        "Notes:\n"
        " • For GTF, exon coordinates are 1-based inclusive; we compute transcript min(start) and max(end) across exons.\n"
        " • TSS/TTS are determined by transcript strand:\n"
        "     + strand: TSS = transcript start, TTS = transcript end\n"
        "     - strand: TSS = transcript end,   TTS = transcript start\n"
        " • Output intervals are 1 bp long by default, but you can expand with --pad."
    )
    parser = argparse.ArgumentParser(description=desc)

    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--gtf',   help='Input GTF with exon features.')
    g.add_argument('--bed12', help='Input BED12 with block structure (exons).')

    parser.add_argument('--site', choices=['tss', 'tts'], required=True,
                        help="Which site to output: 'tss' or 'tts'.")

    parser.add_argument('-o', '--out', required=True,
                        help='Output BED6 file path.')

    parser.add_argument('--name-field',
                        choices=['auto', 'transcript_id', 'gene_id', 'both', 'bed_name'],
                        default='auto',
                        help="How to set the BED name field. Default: auto "
                             "(GTF: transcript_id > gene_id; BED: name column).")

    parser.add_argument('--score', type=int, default=0,
                        help='BED score (0–1000). Default: 0.')

    parser.add_argument('--skip-unknown-strand', action='store_true',
                        help="Skip entries with unknown/ambiguous strand ('.').")

    parser.add_argument('--pad', type=int, default=0,
                        help="Expand sites by this many bases on each side. "
                             "E.g. --pad 50 gives a 101 bp window centered on the site. Default: 0 (1 bp).")

    # Show help if no args when running from CLI; empty args in notebooks.
    if is_interactive():
        args = parser.parse_args('')
    else:
        args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args


def parse_gtf_transcript_bounds(path):
    """
    From a GTF, collect per-transcript bounds across exons.
    Returns dict keyed by (chrom, txid) -> dict with:
      {'strand': '+/-/.', 'min_start_1b': int, 'max_end_1b': int, 'gene_id': str or None}
    """
    tx = {}
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
            m_tx = re.search(r'transcript_id\s+"([^"]+)"', attrs) or re.search(r"transcript_id\s+'([^']+)'", attrs)
            m_gene = re.search(r'gene_id\s+"([^"]+)"', attrs) or re.search(r"gene_id\s+'([^']+)'", attrs)
            if not m_tx:
                if not m_gene:
                    continue
                txid = f"{m_gene.group(1)}:exonset"
            else:
                txid = m_tx.group(1)

            s = int(start)
            e = int(end)
            key = (chrom, txid)
            if key not in tx:
                tx[key] = {
                    'strand': strand,
                    'min_start_1b': s,
                    'max_end_1b': e,
                    'gene_id': (m_gene.group(1) if m_gene else None)
                }
            else:
                rec = tx[key]
                if s < rec['min_start_1b']:
                    rec['min_start_1b'] = s
                if e > rec['max_end_1b']:
                    rec['max_end_1b'] = e
                if rec['gene_id'] is None and m_gene:
                    rec['gene_id'] = m_gene.group(1)
    return tx


def parse_bed12_transcript_bounds(path):
    """
    From a BED12, produce per-line transcript bounds.
    Returns list of dicts with:
      {'chrom','strand','tx_start_0b','tx_end_0b','name'}
    """
    records = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            chrom = parts[0]
            chrom_start = int(parts[1])
            chrom_end   = int(parts[2])
            name  = parts[3] if len(parts) > 3 else '.'
            strand = parts[5] if len(parts) > 5 and parts[5] in ('+','-','.') else '.'
            records.append({
                'chrom': chrom,
                'strand': strand,
                'tx_start_0b': chrom_start,
                'tx_end_0b': chrom_end,
                'name': name
            })
    return records


def choose_name(name_field_mode, chrom, strand, txid=None, gene_id=None, bed_name=None):
    """Resolve BED6 'name' field according to user preference."""
    if name_field_mode == 'bed_name':
        return bed_name if bed_name else '.'
    if name_field_mode == 'transcript_id':
        return txid if txid else (gene_id if gene_id else '.')
    if name_field_mode == 'gene_id':
        return gene_id if gene_id else (txid if txid else '.')
    if name_field_mode == 'both':
        if txid and gene_id:
            return f"{txid}|{gene_id}"
        return txid or gene_id or '.'
    # auto
    if txid or gene_id:
        return txid if txid else gene_id
    return bed_name if bed_name else '.'


def write_bed6(path_out, rows):
    with open(path_out, 'w') as out:
        for r in rows:
            out.write('\t'.join(map(str, r)) + '\n')


########
# MAIN #
########

def main():
    args = parse_args()
    out_path = Path(args.out)

    strand_counts = Counter()
    emitted = 0
    rows = []

    if args.gtf:
        in_path = Path(args.gtf)
        if not in_path.is_file():
            logging.error(f"File not found: {in_path}")
            sys.exit(1)
        logging.info(f"Reading GTF: {in_path}")
        tx_bounds = parse_gtf_transcript_bounds(str(in_path))
        logging.info(f"Collected bounds for {len(tx_bounds)} transcripts")

        for (chrom, txid), rec in tx_bounds.items():
            strand = rec['strand'] if rec['strand'] in ('+','-','.') else '.'
            if args.skip_unknown_strand and strand == '.':
                continue

            tx_start_1b = rec['min_start_1b']
            tx_end_1b   = rec['max_end_1b']

            if args.site == 'tss':
                pos_1b = tx_start_1b if strand == '+' else (tx_end_1b if strand == '-' else tx_start_1b)
            else:  # tts
                pos_1b = tx_end_1b   if strand == '+' else (tx_start_1b if strand == '-' else tx_end_1b)

            center0 = pos_1b - 1
            start0 = max(0, center0 - args.pad)
            end0   = center0 + args.pad + 1

            name = choose_name(args.name_field, chrom, strand, txid=txid, gene_id=rec['gene_id'])
            rows.append((chrom, start0, end0, name, max(0, min(args.score, 1000)), strand))

            emitted += 1
            strand_counts[strand] += 1

    else:
        in_path = Path(args.bed12)
        if not in_path.is_file():
            logging.error(f"File not found: {in_path}")
            sys.exit(1)
        logging.info(f"Reading BED12: {in_path}")
        records = parse_bed12_transcript_bounds(str(in_path))
        logging.info(f"Collected {len(records)} transcript entries")

        for rec in records:
            chrom  = rec['chrom']
            strand = rec['strand'] if rec['strand'] in ('+','-','.') else '.'
            if args.skip_unknown_strand and strand == '.':
                continue

            tx_start_0b = rec['tx_start_0b']
            tx_end_0b   = rec['tx_end_0b']

            if args.site == 'tss':
                pos0 = tx_start_0b if strand == '+' else (tx_end_0b - 1 if strand == '-' else tx_start_0b)
            else:  # tts
                pos0 = tx_end_0b - 1 if strand == '+' else (tx_start_0b if strand == '-' else tx_end_0b - 1)

            center0 = pos0
            start0 = max(0, center0 - args.pad)
            end0   = center0 + args.pad + 1

            name = choose_name(args.name_field, chrom, strand, bed_name=rec['name'])
            rows.append((chrom, start0, end0, name, max(0, min(args.score, 1000)), strand))

            emitted += 1
            strand_counts[strand] += 1

    write_bed6(str(out_path), rows)

    plus = strand_counts.get('+', 0)
    minus = strand_counts.get('-', 0)
    dot = strand_counts.get('.', 0)
    logging.info(f"Wrote {emitted} {args.site.upper()} site(s) to {out_path}")
    logging.info(f"Strand distribution: +={plus}, -={minus}, .={dot}")

if __name__ == '__main__':
    main()
