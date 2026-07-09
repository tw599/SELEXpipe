#!/usr/bin/env python3
"""Trim filtered FASTA sequences to fixed-length k-mers (default: 40-mers).

By default each retained sequence is written as a plain one-sequence-per-line
``.seq`` file. Use ``--fasta`` to emit FASTA instead. The extraction window is
controlled by ``--start`` / ``--center``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Trim FASTA sequences to fixed-length k-mers (e.g. 40-mers)."
    )
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file (typically length-filtered)")
    parser.add_argument("output_path", type=Path, help="Output .seq (or FASTA) path")
    parser.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        default=40,
        help="K-mer length to extract (default: 40)",
    )
    start_group = parser.add_mutually_exclusive_group()
    start_group.add_argument(
        "-s",
        "--start",
        type=int,
        default=0,
        help="0-based start offset for the k-mer window (default: 0)",
    )
    start_group.add_argument(
        "--center",
        action="store_true",
        help="Center the k-mer window within each sequence",
    )
    parser.add_argument(
        "--fasta",
        action="store_true",
        help="Write FASTA output instead of one sequence per line",
    )
    return parser.parse_args()


def iter_fasta(path: Path):
    header: str | None = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).replace(" ", "").replace("\t", "")
                header = line[1:].strip() or line
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"Malformed FASTA in {path}: sequence before header")
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks).replace(" ", "").replace("\t", "")


def trim_to_kmers(
    input_fasta: Path,
    output_path: Path,
    kmer_size: int,
    start: int,
    center: bool,
    as_fasta: bool,
) -> tuple[int, int]:
    if kmer_size <= 0:
        raise ValueError(f"k-mer size must be positive, got {kmer_size}")
    if start < 0:
        raise ValueError(f"start offset must be >= 0, got {start}")

    total = 0
    written = 0

    with output_path.open("w", encoding="utf-8") as out:
        for header, sequence in iter_fasta(input_fasta):
            total += 1
            if len(sequence) < kmer_size:
                continue

            if center:
                window_start = (len(sequence) - kmer_size) // 2
            else:
                window_start = start

            window_end = window_start + kmer_size
            if window_start < 0 or window_end > len(sequence):
                continue

            kmer = sequence[window_start:window_end]
            if as_fasta:
                out.write(f">{header}\n")
                for i in range(0, len(kmer), 80):
                    out.write(kmer[i : i + 80] + "\n")
            else:
                out.write(kmer + "\n")
            written += 1

    return total, written


def main() -> int:
    args = parse_args()
    if not args.input_fasta.is_file():
        print(f"error: input FASTA not found: {args.input_fasta}", file=sys.stderr)
        return 1

    try:
        total, written = trim_to_kmers(
            args.input_fasta,
            args.output_path,
            args.kmer_size,
            args.start,
            args.center,
            args.fasta,
        )
    except Exception as exc:  # noqa: BLE001 - surface parse/IO errors to CLI
        print(f"error: {exc}", file=sys.stderr)
        return 1

    mode = "centered" if args.center else f"start={args.start}"
    print(
        f"trim_to_kmer: wrote {written}/{total} {args.kmer_size}-mers ({mode}) "
        f"-> {args.output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
