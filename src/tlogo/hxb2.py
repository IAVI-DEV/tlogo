"""HxB2 coordinate mapping for HIV Env gp160 protein alignments.

Maps alignment columns to HxB2 amino acid positions and Env structural
regions.  HxB2 is the standard HIV-1 reference numbering -- since Env
varies wildly in V-loop length, HxB2 positions provide a universal
coordinate system for cross-strain comparisons.

    >>> get_env_region(1)
    'SP'
    >>> get_env_region(131)
    'V1'
    >>> get_env_region(296)
    'V3'
    >>> get_env_region(512)
    'gp41'
    >>> get_env_region(0) is None
    True
    >>> get_env_region(857) is None
    True
"""

from __future__ import annotations

from dataclasses import dataclass

from tlogo.constants import ENV_REGIONS

# Pre-build lookup for fast position -> region resolution
_REGION_LOOKUP: dict[int, str] = {}
for _name, _start, _end in ENV_REGIONS:
    for _pos in range(_start, _end + 1):
        _REGION_LOOKUP[_pos] = _name


def get_env_region(hxb2_aa_pos: int) -> str | None:
    """Return the Env structural region name for an HxB2 AA position.

    Args:
        hxb2_aa_pos: 1-based HxB2 amino acid position.

    Returns:
        Region name (e.g. "V3", "C1", "gp41") or None if out of range.

    >>> get_env_region(150)
    'V1'
    >>> get_env_region(300)
    'V3'
    >>> get_env_region(500)
    'C5'
    """
    return _REGION_LOOKUP.get(hxb2_aa_pos)


@dataclass
class HxB2Position:
    """One alignment column mapped to HxB2 coordinates.

    Attributes:
        alignment_col: 0-based column index in the alignment.
        hxb2_aa_pos: 1-based HxB2 AA position, or None if this column
            is an insertion relative to HxB2 (i.e. HxB2 has a gap here).
        region: Env region name, or None for insertion columns.
        hxb2_residue: The HxB2 residue at this position ("-" for gaps).
    """

    alignment_col: int
    hxb2_aa_pos: int | None
    region: str | None
    hxb2_residue: str


def build_hxb2_map(
    alignment: dict[str, str],
    hxb2_id: str = "HxB2",
) -> list[HxB2Position]:
    """Walk the HxB2 row of a protein alignment and map every column.

    Non-gap positions in the HxB2 sequence are numbered sequentially
    (1-based).  Gap positions (insertions relative to HxB2) get
    ``hxb2_aa_pos=None``.

    Args:
        alignment: {seq_id: aligned_sequence} dict.
        hxb2_id: Sequence ID of the HxB2 reference in the alignment.

    Returns:
        List of HxB2Position, one per alignment column.

    Raises:
        ValueError: If hxb2_id is not found in the alignment.
    """
    hxb2_seq = alignment.get(hxb2_id)
    if hxb2_seq is None:
        raise ValueError(
            f"HxB2 sequence '{hxb2_id}' not found in alignment. "
            f"Available IDs: {list(alignment.keys())[:10]}"
        )

    positions: list[HxB2Position] = []
    aa_counter = 0

    for col_idx, residue in enumerate(hxb2_seq):
        if residue == "-":
            positions.append(HxB2Position(
                alignment_col=col_idx,
                hxb2_aa_pos=None,
                region=None,
                hxb2_residue="-",
            ))
        else:
            aa_counter += 1
            positions.append(HxB2Position(
                alignment_col=col_idx,
                hxb2_aa_pos=aa_counter,
                region=get_env_region(aa_counter),
                hxb2_residue=residue,
            ))

    return positions


def build_reverse_hxb2_map(
    hxb2_map: list[HxB2Position],
) -> dict[int, int]:
    """Build mapping from HxB2 AA position -> alignment column index.

    Only includes non-None hxb2_aa_pos entries.

    >>> pos = [HxB2Position(0, 10, "C1", "A"), HxB2Position(1, None, None, "-"), HxB2Position(2, 11, "C1", "B")]
    >>> build_reverse_hxb2_map(pos)
    {10: 0, 11: 2}
    """
    return {
        p.hxb2_aa_pos: p.alignment_col
        for p in hxb2_map
        if p.hxb2_aa_pos is not None
    }
