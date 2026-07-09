#!/usr/bin/env python3
"""Filter a FASTA file to sequences of an expected length.

Retains only records whose sequence length equals ``--length`` (default: 101),
matching the expected merged-read size after adapter trimming and read merging.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Keep FASTA records whose sequence length equals the expected read size."
    )
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_fasta", type=Path, help="Output FASTA file (filtered)")
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        default=101,
        help="Expected sequence length to retain (default: 101)",
    )
    return parser.parse_args()


def filter_fasta(input_fasta: Path, output_fasta: Path, expected_length: int) -> tuple[int, int]:
    if expected_length <= 0:
        raise ValueError(f"expected length must be positive, got {expected_length}")

    kept = 0
    total = 0
    header: str | None = None
    seq_chunks: list[str] = []

    def flush() -> None:
        nonlocal kept, total, header, seq_chunks
        if header is None:
            return
        sequence = "".join(seq_chunks).replace(" ", "").replace("\t", "")
        total += 1
        if len(sequence) == expected_length:
            output.write(f"{header}\n")
            # wrap at 80 columns for readability
            for i in range(0, len(sequence), 80):
                output.write(sequence[i : i + 80] + "\n")
            kept += 1
        header = None
        seq_chunks = []

    with input_fasta.open("r", encoding="utf-8") as handle, output_fasta.open(
        "w", encoding="utf-8"
    ) as output:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"Malformed FASTA in {input_fasta}: sequence before header")
                seq_chunks.append(line.strip())
        flush()

    return total, kept


def main() -> int:
    args = parse_args()
    if not args.input_fasta.is_file():
        print(f"error: input FASTA not found: {args.input_fasta}", file=sys.stderr)
        return 1

    try:
        total, kept = filter_fasta(args.input_fasta, args.output_fasta, args.length)
    except Exception as exc:  # noqa: BLE001 - surface parse/IO errors to CLI
        print(f"error: {exc}", file=sys.stderr)
        return 1

    discarded = total - kept
    print(
        f"filter_by_length: kept {kept}/{total} sequences of length {args.length} "
        f"(discarded {discarded}) -> {args.output_fasta}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
