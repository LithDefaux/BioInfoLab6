#!/usr/bin/env python3
"""
Bioinformatics lab: simulating and detecting transposable elements (TEs).

Part 1 - build an artificial DNA sequence (200-400 bases) with 3-4 TEs.
Part 2 - detect the (start, end) of every TE in that sequence.

Model used for the TEs ("cut-and-paste" DNA transposons):

    TSD | TIR ........ internal region ........ TIR' | TSD
        <------------------- TE ------------------->

  * TIR  = Terminal Inverted Repeat. The right end of the element is the
           reverse complement of the left end. This is the "signature"
           that lets us find a TE without knowing its sequence.
  * TSD  = Target Site Duplication. When a TE inserts, the host DNA at the
           insertion point (here 5 bases) is duplicated on both sides.

Two detectors are implemented:
  A) de novo detector  - looks for inverted repeats (TIRs) and confirms them
                         with the direct repeats (TSDs) on the flanks. Does
                         NOT need to know the TE sequences.
  B) library detector  - plain pattern search of known TE sequences
                         (both strands). Works as a sanity check.

Usage:
    python te_lab.py                 # random run
    python te_lab.py --seed 42       # reproducible run
    python te_lab.py --fasta out.fa  # also save the sequence as FASTA
"""

import argparse
import random

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

TSD_LEN = 5          # length of the target site duplication (5 bases, as on the board)
MIN_TIR, MAX_TIR = 10, 12
MIN_TE_LEN, MAX_TE_LEN = 36, 50


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def random_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


# --------------------------------------------------------------------------
# PART 1 - generate the artificial sequence
# --------------------------------------------------------------------------
def make_te(rng):
    """Create one TE: TIR + random internal region + reverse complement of TIR."""
    tir_len = rng.randint(MIN_TIR, MAX_TIR)
    total = rng.randint(MIN_TE_LEN, MAX_TE_LEN)
    tir = random_dna(rng, tir_len)
    internal = random_dna(rng, total - 2 * tir_len)
    return tir + internal + revcomp(tir)


def generate_sequence(rng, n_te=None):
    """
    Return (sequence, truth) where truth is a list of dicts with the real
    position of every inserted TE (0-based, end inclusive).
    """
    n_te = n_te or rng.randint(3, 4)
    background = random_dna(rng, rng.randint(110, 160))
    tes = [make_te(rng) for _ in range(n_te)]

    # choose insertion points on the background: spaced out so TEs never
    # overlap or nest, and with a TSD whose edges would not accidentally
    # extend the TIR (keeps the TE boundaries unambiguous).
    while True:
        points = sorted(rng.sample(range(10, len(background) - 10), n_te))
        spaced = all(b - a >= TSD_LEN + 8 for a, b in zip(points, points[1:]))
        clean = all(
            background[p + TSD_LEN - 1] != COMP[background[p]] for p in points
        )
        if spaced and clean:
            break

    pieces, truth, cursor, length = [], [], 0, 0
    for k, (p, te) in enumerate(zip(points, tes), start=1):
        tsd = background[p:p + TSD_LEN]

        chunk = background[cursor:p] + tsd          # host DNA + left TSD
        pieces.append(chunk)
        length += len(chunk)

        truth.append({"name": f"TE{k}", "start": length,
                      "end": length + len(te) - 1, "seq": te, "tsd": tsd})
        pieces.append(te)
        length += len(te)

        pieces.append(tsd)                          # right TSD
        length += TSD_LEN
        cursor = p + TSD_LEN

    pieces.append(background[cursor:])
    seq = "".join(pieces)
    assert 200 <= len(seq) <= 400, len(seq)
    return seq, truth


# --------------------------------------------------------------------------
# PART 2A - de novo detector (terminal inverted repeats)
# --------------------------------------------------------------------------
def find_tir_candidates(seq, min_tir=MIN_TIR, min_len=30, max_len=60):
    """
    For every window [i, j] check whether its two ends are inverted repeats:
        seq[i + k] == complement(seq[j - k])     for k = 0 .. L-1
    Keep (i, j, L) if L >= min_tir and the match cannot be extended outward
    (so we report the outermost, i.e. true, boundaries of the element).
    """
    n = len(seq)
    hits = []
    for i in range(n - min_len + 1):
        for j in range(i + min_len - 1, min(i + max_len, n)):
            # maximality: skip if the pair just outside also matches
            if i > 0 and j < n - 1 and seq[i - 1] == COMP[seq[j + 1]]:
                continue
            L = 0
            while i + L < j - L and seq[i + L] == COMP[seq[j - L]]:
                L += 1
            if L >= min_tir:
                hits.append((i, j, L))
    return hits


def flanking_tsd(seq, i, j, k=TSD_LEN):
    """Return the direct repeat (TSD) around the element [i, j], or None."""
    if i < k or j + k >= len(seq):
        return None
    left, right = seq[i - k:i], seq[j + 1:j + 1 + k]
    return left if left == right else None


def detect_denovo(seq, require_tsd=True):
    """
    Candidates = inverted repeats (TIRs). With require_tsd=True a candidate
    is accepted only if it is also flanked by direct repeats (TSDs), exactly
    the two signatures drawn on the board. Overlaps: longest TIR wins.
    Returns a sorted list of (start, end, tir_len, tsd).
    """
    candidates = sorted(find_tir_candidates(seq),
                        key=lambda h: (-h[2], -(h[1] - h[0])))
    chosen = []
    for i, j, L in candidates:
        tsd = flanking_tsd(seq, i, j)
        if require_tsd and tsd is None:
            continue
        if all(j < a or i > b for a, b, _, _ in chosen):
            chosen.append((i, j, L, tsd))
    return sorted(chosen)


# --------------------------------------------------------------------------
# PART 2B - library detector (known TE sequences)
# --------------------------------------------------------------------------
def detect_with_library(seq, library):
    """library: {name: sequence}. Searches forward and reverse strand."""
    found = []
    for name, te in library.items():
        for strand, query in (("+", te), ("-", revcomp(te))):
            pos = seq.find(query)
            while pos != -1:
                found.append((name, strand, pos, pos + len(te) - 1))
                pos = seq.find(query, pos + 1)
    return sorted(found, key=lambda f: f[2])


# --------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------
def annotate(seq, truth):
    """Return the sequence with TE regions in lowercase."""
    s = list(seq)
    for t in truth:
        for x in range(t["start"], t["end"] + 1):
            s[x] = s[x].lower()
    return "".join(s)


def print_wrapped(s, width=60):
    for i in range(0, len(s), width):
        print(f"{i + 1:>4}  {s[i:i + width]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--fasta", help="save the sequence to this FASTA file")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    # ---- Part 1
    seq, truth = generate_sequence(rng)
    print(f"Artificial sequence: {len(seq)} bases, {len(truth)} transposable elements")
    print("(lowercase = inserted TE)\n")
    print_wrapped(annotate(seq, truth))

    if args.fasta:
        with open(args.fasta, "w") as f:
            f.write(">artificial_dna_with_TEs\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
        print(f"\nSaved to {args.fasta}")

    print("\nGround truth (1-based, inclusive):")
    for t in truth:
        print(f"  {t['name']}: start={t['start'] + 1:>3}  end={t['end'] + 1:>3}  "
              f"len={len(t['seq'])}  TSD={t['tsd']}")

    # ---- Part 2A
    hits = detect_denovo(seq)
    print("\nDetector A - de novo (inverted repeats + direct repeats):")
    for n, (i, j, L, tsd) in enumerate(hits, start=1):
        print(f"  candidate {n}: start={i + 1:>3}  end={j + 1:>3}  "
              f"len={j - i + 1}  TIR>={L}  TSD={tsd}")

    real = {(t["start"], t["end"]) for t in truth}
    found = {(i, j) for i, j, _, _ in hits}
    print(f"  -> exact matches: {len(real & found)}/{len(real)} | "
          f"false positives: {len(found - real)}")

    # ---- Part 2B
    library = {t["name"]: t["seq"] for t in truth}
    print("\nDetector B - known-library search (both strands):")
    for name, strand, s, e in detect_with_library(seq, library):
        print(f"  {name} ({strand}): start={s + 1:>3}  end={e + 1:>3}")


if __name__ == "__main__":
    main()