#!/usr/bin/env python3
"""Convert FASTA reads into a SEQ file of fixed-length lines (default: 40-mers).

Each FASTA record is split into consecutive chunks of ``--length`` bases.
Trailing chunks shorter than ``--length`` are padded on the right with random
A/T/C/G bases. Output is one sequence per line (no FASTA headers).
"""

from __future__ import annotations

import argparse
import random
import sys


def pad_sequence(seq: str, length: int = 40) -> str:
    if len(seq) < length:
        padding = "".join(random.choices("ATCG", k=length - len(seq)))
        return seq + padding
    return seq


def process_fasta_to_seq(
    input_fasta: str,
    output_seq: str,
    line_length: int = 40,
) -> dict[str, int]:
    if line_length <= 0:
        raise ValueError(f"length must be positive, got {line_length}")

    n_records = 0
    n_lines = 0
    n_padded = 0

    with open(input_fasta, "r") as infile, open(output_seq, "w") as outfile:
        current_seq = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    n_records += 1
                    for i in range(0, len(current_seq), line_length):
                        chunk = current_seq[i : i + line_length]
                        if len(chunk) < line_length:
                            n_padded += 1
                        outfile.write(pad_sequence(chunk, line_length) + "\n")
                        n_lines += 1
                    current_seq = ""
            else:
                current_seq += line.upper()

        # Process the last read
        if current_seq:
            n_records += 1
            for i in range(0, len(current_seq), line_length):
                chunk = current_seq[i : i + line_length]
                if len(chunk) < line_length:
                    n_padded += 1
                outfile.write(pad_sequence(chunk, line_length) + "\n")
                n_lines += 1

    return {
        "records": n_records,
        "lines": n_lines,
        "padded_chunks": n_padded,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Convert FASTA reads into a SEQ file with fixed-length lines "
            "(default 40 bases), padded with random bases if needed."
        )
    )
    parser.add_argument("input", help="Path to the input FASTA file")
    parser.add_argument("output", help="Path to the output SEQ file")
    parser.add_argument(
        "--length",
        type=int,
        default=40,
        help="Line length for each sequence fragment (default: 40)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility (optional)",
    )

    args = parser.parse_args(argv)

    if args.seed is not None:
        random.seed(args.seed)

    try:
        stats = process_fasta_to_seq(args.input, args.output, args.length)
    except Exception as exc:  # noqa: BLE001 - surface parse/IO errors to CLI
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(
        f"[OK] FASTA records:   {stats['records']}\n"
        f"[OK] SEQ lines:       {stats['lines']} "
        f"(length={args.length})\n"
        f"[OK] Padded chunks:   {stats['padded_chunks']}\n"
        f"[OK] Wrote:           {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
