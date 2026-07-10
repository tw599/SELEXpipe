#!/usr/bin/env python3
"""Normalize FASTA reads to a user-specified fixed length.

Keeps only reads within an allowed input-length window, then pads shorter
reads (random A/T/C/G) or trims longer reads so every retained sequence is
exactly ``--target_length`` bases.

Expected merged-read lengths for common ligand designs:

| Ligand design | Expected length |
|---------------|-----------------|
| lig147        | 101 bp          |
| lig200        | 154 bp          |
| ligN40        | 187 bp          |
| ligN70        | 247 bp          |

Use ``--ligand`` to select a design, or ``--target_length`` for an explicit value.
"""

from __future__ import annotations

import argparse
import random
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Expected merged-read lengths for supported ligand designs.
LIGAND_LENGTHS: dict[str, int] = {
    "lig147": 101,
    "lig200": 154,
    "ligN40": 187,
    "ligN70": 247,
}


def normalize_fasta_to_fixed_length(
    input_fasta: str,
    output_fasta: str,
    target_length: int,
    min_input_len: int,
    max_input_len: int,
    trim_from: str = "right",
) -> dict[str, int]:
    """
    Keep only reads within an allowed input-length window, then normalize all
    kept reads to exactly ``target_length`` by padding shorter reads or
    trimming longer reads.

    Args:
        input_fasta: Input FASTA file
        output_fasta: Output FASTA file
        target_length: Final required length for all retained reads
        min_input_len: Minimum input length to retain
        max_input_len: Maximum input length to retain
        trim_from: ``right``, ``left``, or ``center``
    """
    if target_length <= 0:
        raise ValueError(f"target_length must be positive, got {target_length}")
    if min_input_len <= 0:
        raise ValueError(f"min_input_len must be positive, got {min_input_len}")
    if max_input_len < min_input_len:
        raise ValueError(
            f"max_input_len ({max_input_len}) must be >= min_input_len ({min_input_len})"
        )
    if not (min_input_len <= target_length <= max_input_len):
        raise ValueError(
            f"target_length ({target_length}) must lie within "
            f"[{min_input_len}, {max_input_len}]"
        )
    if trim_from not in {"right", "left", "center"}:
        raise ValueError("--trim_from must be one of: right, left, center")

    bases = ["A", "T", "C", "G"]

    n_total = 0
    n_kept = 0
    n_padded = 0
    n_trimmed = 0
    n_discarded = 0

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            n_total += 1
            sequence = str(record.seq).upper()
            length = len(sequence)

            # Discard reads outside allowed input range
            if length < min_input_len or length > max_input_len:
                n_discarded += 1
                continue

            # Pad shorter reads
            if length < target_length:
                padding_length = target_length - length
                sequence = sequence + "".join(random.choices(bases, k=padding_length))
                n_padded += 1

            # Trim longer reads
            elif length > target_length:
                if trim_from == "right":
                    sequence = sequence[:target_length]
                elif trim_from == "left":
                    sequence = sequence[-target_length:]
                else:  # center
                    extra = length - target_length
                    left_trim = extra // 2
                    sequence = sequence[left_trim : left_trim + target_length]
                n_trimmed += 1

            # Keep exact-length reads unchanged
            assert len(sequence) == target_length

            adjusted_record = SeqRecord(
                Seq(sequence),
                id=record.id,
                description=record.description,
            )
            SeqIO.write(adjusted_record, outfile, "fasta")
            n_kept += 1

    return {
        "total": n_total,
        "kept": n_kept,
        "padded": n_padded,
        "trimmed": n_trimmed,
        "discarded": n_discarded,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ligand_help = ", ".join(
        f"{name}={length} bp" for name, length in LIGAND_LENGTHS.items()
    )
    parser = argparse.ArgumentParser(
        description=(
            "Filter FASTA reads by an input-length window, then normalize all "
            "retained reads to a user-specified fixed length (pad or trim)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Ligand design expected lengths:\n"
            "  lig147  101 bp\n"
            "  lig200  154 bp\n"
            "  ligN40  187 bp\n"
            "  ligN70  247 bp\n"
        ),
    )
    parser.add_argument("input_fasta", type=str, help="Input FASTA file")
    parser.add_argument("output_fasta", type=str, help="Output FASTA file")
    length_group = parser.add_mutually_exclusive_group(required=True)
    length_group.add_argument(
        "--target_length",
        type=int,
        default=None,
        help="Final required read length (no hard-coded default)",
    )
    length_group.add_argument(
        "--ligand",
        choices=sorted(LIGAND_LENGTHS),
        default=None,
        help=f"Ligand design key used to set target length ({ligand_help})",
    )
    parser.add_argument(
        "--min_input_len",
        type=int,
        default=None,
        help=(
            "Minimum input read length to retain "
            "(default: same as --target_length)"
        ),
    )
    parser.add_argument(
        "--max_input_len",
        type=int,
        default=None,
        help=(
            "Maximum input read length to retain "
            "(default: same as --target_length)"
        ),
    )
    parser.add_argument(
        "--trim_from",
        choices=["right", "left", "center"],
        default="right",
        help="How to trim reads longer than target_length (default: right)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible padding",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.ligand is not None:
        target_length = LIGAND_LENGTHS[args.ligand]
        ligand_label = args.ligand
    else:
        target_length = args.target_length
        ligand_label = None

    min_input_len = (
        args.min_input_len if args.min_input_len is not None else target_length
    )
    max_input_len = (
        args.max_input_len if args.max_input_len is not None else target_length
    )

    if args.seed is not None:
        random.seed(args.seed)

    try:
        stats = normalize_fasta_to_fixed_length(
            input_fasta=args.input_fasta,
            output_fasta=args.output_fasta,
            target_length=target_length,
            min_input_len=min_input_len,
            max_input_len=max_input_len,
            trim_from=args.trim_from,
        )
    except Exception as exc:  # noqa: BLE001 - surface errors to CLI
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if ligand_label is not None:
        print(f"[OK] Ligand design:  {ligand_label}")
    print(f"[OK] Target length:  {target_length}")
    print(f"[OK] Input window:   {min_input_len}-{max_input_len}")
    print(f"[OK] Total reads:    {stats['total']}")
    print(f"[OK] Reads kept:     {stats['kept']}")
    print(f"[OK] Reads padded:   {stats['padded']}")
    print(f"[OK] Reads trimmed:  {stats['trimmed']}")
    print(f"[OK] Reads dropped:  {stats['discarded']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
