# Transposable Element Simulation and Detection

A bioinformatics lab project in Python. It builds an artificial DNA sequence
containing simulated transposable elements (TEs) and then detects the start and
end position of each element.

## Task

1. Make an artificial DNA sequence of 200 to 400 bases in which 3 to 4
   transposable elements are simulated.
2. Implement a software application that detects the positions (start, end) of
   these elements within the created sequence.

## How it works

### Simulation

Each element is modelled as a cut-and-paste DNA transposon:

```
TSD | TIR ........ internal region ........ TIR' | TSD
    <------------------- TE ------------------->
```

- **TIR (terminal inverted repeat):** the right end of the element is the
  reverse complement of its left end.
- **TSD (target site duplication):** the 5-base target site in the host DNA is
  duplicated on both sides of the inserted element (the "direct repeats").

The generator creates a random background sequence, inserts 3 or 4 distinct
elements at non-overlapping positions, and keeps the true positions as ground
truth. The final sequence length is always between 200 and 400 bases.

### Detection

**Detector A (de novo).** It does not need to know the element sequences. It
scans every window of 30 to 60 bases and keeps those whose two ends are inverted
repeats (at least 10 bases). A candidate is accepted only if it is also flanked
by matching direct repeats (TSDs). The match is extended outward as far as
possible so the reported boundaries are the outermost ones. When candidates
overlap, the one with the longest TIR is kept.

**Detector B (library search).** A plain search of known element sequences on
both strands. It serves as a sanity check against Detector A.

The program compares Detector A with the ground truth and prints the number of
exact matches and false positives.

## Requirements

- Python 3.8 or newer
- No third-party packages

## Usage

```bash
python te_lab.py                  # random run
python te_lab.py --seed 42        # reproducible run
python te_lab.py --fasta out.fa   # also save the sequence in FASTA format
```

Positions are printed 1-based and inclusive.

## Example output

```
Artificial sequence: 303 bases, 3 transposable elements

Ground truth (1-based, inclusive):
  TE1: start= 29  end= 75  len=47  TSD=TCTAC
  TE2: start=123  end=168  len=46  TSD=TAGTC
  TE3: start=235  end=283  len=49  TSD=TAAGT

Detector A - de novo (inverted repeats + direct repeats):
  candidate 1: start= 29  end= 75  len=47  TIR>=11  TSD=TCTAC
  candidate 2: start=123  end=168  len=46  TIR>=10  TSD=TAGTC
  candidate 3: start=235  end=283  len=49  TIR>=11  TSD=TAAGT
  -> exact matches: 3/3 | false positives: 0
```

(Run with `--seed 1`; the sequence itself is also printed, with the inserted
elements in lowercase.)

## Evaluation

Over 1,000 random runs (3,497 elements) Detector A found every element at the
exact start and end positions, with no false positives.

## Limitations

- Detector A needs perfect inverted repeats. Mutated copies would require
  mismatch-tolerant matching.
- Elements are inserted once each and never nest inside each other.
- The sequence is single-stranded; the staggered cut and gap filling that
  produce the target site duplication are not simulated step by step.

## Project structure

```
.
├── te_lab.py     # generator and detectors
├── README.md
└── .gitignore
```
