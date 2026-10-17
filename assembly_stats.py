#!/usr/bin/env python3
"""
assembly_stats.py - summarise de novo assembly FASTA files (Velvet contigs.fa / Oases transcripts.fa)

Usage:
  assembly_stats.py header
      -> prints the TSV header line for the summary table
  assembly_stats.py summary <tool> <k> <fasta> [wall_seconds] [max_rss_kb]
      -> prints one TSV row: tool, k, #sequences, total bp, longest, N50, L50,
         #sequences >= 1 kb, #loci (Oases only), wall time (s), peak memory (MB)
  assembly_stats.py lengths <fasta>
      -> prints every sequence name with its length (bp), longest first

Lengths are computed from the sequence itself, NOT from the FASTA header.
(Velvet's "NODE_x_length_L" header reports L in k-mers; true length = L + k - 1.)
"""
import sys


def read_lengths(path):
    """Return a list of (name, length) for every record in a FASTA file."""
    records = []
    name, n = None, 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, n))
                name, n = line[1:].split()[0], 0
            else:
                n += len(line)
    if name is not None:
        records.append((name, n))
    return records


def n50_l50(lengths):
    """N50 = length L such that sequences >= L cover half the assembly; L50 = how many of them."""
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    running = 0
    for i, length in enumerate(s, 1):
        running += length
        if running >= half:
            return length, i
    return 0, 0


HEADER = ["tool", "k", "num_seqs", "total_bp", "longest_bp", "N50", "L50",
          "num_ge_1kb", "num_loci", "wall_s", "peak_mem_MB"]


def summary(tool, k, fasta, wall="NA", rss_kb="NA"):
    recs = read_lengths(fasta)
    lens = [l for _, l in recs]
    n50, l50 = n50_l50(lens)
    if any(name.startswith("Locus_") for name, _ in recs):
        loci = len({name.split("_Transcript")[0] for name, _ in recs})
    else:
        loci = "NA"
    mem = f"{int(rss_kb) / 1024:.0f}" if str(rss_kb).isdigit() else "NA"
    row = [tool, k, len(lens), sum(lens), max(lens, default=0), n50, l50,
           sum(1 for l in lens if l >= 1000), loci, wall, mem]
    print("\t".join(str(x) for x in row))


def lengths(fasta):
    print("name\tlength_bp")
    for name, l in sorted(read_lengths(fasta), key=lambda r: -r[1]):
        print(f"{name}\t{l}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "header":
        print("\t".join(HEADER))
    elif cmd == "summary" and len(sys.argv) >= 5:
        summary(*sys.argv[2:7])
    elif cmd == "lengths" and len(sys.argv) == 3:
        lengths(sys.argv[2])
    else:
        sys.exit(__doc__)
