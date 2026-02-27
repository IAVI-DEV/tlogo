"""FASTA I/O and animal name extraction from sequence IDs."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def parse_alignment(fasta_path: Path) -> dict[str, str]:
    """Parse aligned FASTA into {seq_id: aligned_sequence} dict.

    Args:
        fasta_path: Path to aligned FASTA file.

    Returns:
        Ordered dict of seq_id -> aligned sequence string.
    """
    alignment: dict[str, str] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        alignment[rec.id] = str(rec.seq)
    return alignment


def extract_animal_name(
    seq_id: str,
    delimiter: str = "_",
    field: int = 0,
) -> str:
    """Extract animal/group name from a sequence ID.

    Splits seq_id on delimiter and returns the specified field.

    Args:
        seq_id: The FASTA sequence identifier.
        delimiter: Character to split on.
        field: 0-based index of the field to extract.

    Returns:
        The extracted animal name, or the full seq_id if splitting fails.

    >>> extract_animal_name("rh2856_wk8_001")
    'rh2856'
    >>> extract_animal_name("rh2856_wk8_001", field=1)
    'wk8'
    >>> extract_animal_name("no-delimiter")
    'no-delimiter'
    """
    parts = seq_id.split(delimiter)
    if field < len(parts):
        return parts[field]
    return seq_id


def group_sequences_by_animal(
    alignment: dict[str, str],
    exclude_ids: set[str],
    delimiter: str = "_",
    field: int = 0,
) -> dict[str, list[str]]:
    """Group alignment sequences by animal name.

    Excludes specified IDs (e.g. HxB2, reference) and groups remaining
    sequences using extract_animal_name.

    Args:
        alignment: {seq_id: aligned_seq} dict.
        exclude_ids: Sequence IDs to skip (references).
        delimiter: Passed to extract_animal_name.
        field: Passed to extract_animal_name.

    Returns:
        {animal_name: [aligned_seq, ...]} dict.
    """
    groups: dict[str, list[str]] = defaultdict(list)
    for seq_id, seq in alignment.items():
        if seq_id in exclude_ids:
            continue
        animal = extract_animal_name(seq_id, delimiter, field)
        groups[animal].append(seq)
    return dict(groups)


def sort_animal_groups(
    animal_names: list[str],
    self_name: str | None = None,
) -> list[str]:
    """Sort animal groups: self first, recombinants, then alphabetical.

    Args:
        animal_names: List of animal name strings.
        self_name: Optional name to sort first (e.g. lineage name).

    Returns:
        Sorted list.

    >>> sort_animal_groups(["rh2856", "CH505", "Rec5050", "r17042"], "CH505")
    ['CH505', 'Rec5050', 'r17042', 'rh2856']
    >>> sort_animal_groups(["beta", "alpha", "RecX"])
    ['RecX', 'alpha', 'beta']
    """
    self_group: list[str] = []
    rec_group: list[str] = []
    other_group: list[str] = []

    for name in animal_names:
        if self_name is not None and name == self_name:
            self_group.append(name)
        elif name.startswith("Rec"):
            rec_group.append(name)
        else:
            other_group.append(name)

    return self_group + sorted(rec_group) + sorted(other_group)


def detect_reference_id(
    alignment: dict[str, str],
    ref_id: str | None = None,
    hxb2_id: str = "HxB2",
) -> str | None:
    """Detect or validate the parental reference sequence ID.

    If ref_id is provided, validates it exists. Otherwise uses heuristics:
    1. Look for a seq_id ending with '_ref'.
    2. Fall back to None (only HxB2 excluded).

    Args:
        alignment: The alignment dict.
        ref_id: Explicit reference ID, or None for auto-detection.
        hxb2_id: HxB2 ID to skip during auto-detection.

    Returns:
        The reference sequence ID, or None if not found.
    """
    if ref_id is not None:
        if ref_id in alignment:
            return ref_id
        return None

    for seq_id in alignment:
        if seq_id == hxb2_id:
            continue
        if seq_id.endswith("_ref"):
            return seq_id

    return None
