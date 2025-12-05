#!/usr/bin/env python3
import sys
import textwrap

if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} ALIGNMENT.aln RECOMB.gff OUTPUT.fasta")

aln_path = sys.argv[1]
gff_path = sys.argv[2]
out_path = sys.argv[3]

# ------------------ Read alignment ------------------
seqs = {}
order = []

with open(aln_path) as f:
    name, buf = None, []
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
                order.append(name)
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())

    if name is not None:
        seqs[name] = "".join(buf)
        order.append(name)

if not seqs:
    sys.exit("ERROR: alignment file is empty or malformed.")

L = len(next(iter(seqs.values())))
if any(len(s) != L for s in seqs.values()):
    sys.exit("ERROR: alignment sequences have unequal lengths.")

# ------------------ Read GFF intervals ------------------
intervals = []

with open(gff_path) as g:
    for line in g:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            continue
        if start > end:
            start, end = end, start
        intervals.append((start, end))

# ------------------ Merge intervals ------------------
def merge_intervals(ranges):
    if not ranges:
        return []
    ranges = sorted(ranges)
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]

intervals = merge_intervals(intervals)

# ------------------ Build mask array ------------------
mask = bytearray(b"\x00") * L

for s, e in intervals:
    s0 = max(0, s - 1)
    e0 = min(L - 1, e - 1)
    if s0 <= e0:
        mask[s0:e0 + 1] = b"\x01" * (e0 - s0 + 1)

# ------------------ Apply masking ------------------
with open(out_path, "w") as out:
    for name in order:
        seq = list(seqs[name])
        for i, m in enumerate(mask):
            if m and seq[i] != "-":
                seq[i] = "N"
        masked = "".join(seq)
        out.write(f">{name}\n")
        for chunk in textwrap.wrap(masked, 60):
            out.write(chunk + "\n")
