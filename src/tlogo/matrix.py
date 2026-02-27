"""Position resolution, gap filtering, and logomaker matrix building.

    >>> parse_positions("130,160,332")
    [130, 160, 332]
    >>> parse_positions(" 50 , 100 , 200 ")
    [50, 100, 200]
    >>> parse_regions("V3")
    ['V3']
    >>> parse_regions("V1,V2,V3")
    ['V1', 'V2', 'V3']
    >>> compute_gap_fraction(["---", "---", "---"], 0)
    1.0
    >>> compute_gap_fraction(["ABC", "DEF", "GHI"], 0)
    0.0
    >>> compute_gap_fraction(["A-C", "-BC", "ABC", "A-C"], 1)
    0.5
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import logomaker
import pandas as pd

from tlogo.constants import ENV_REGIONS, MAX_GAP_FRACTION
from tlogo.hxb2 import HxB2Position, get_env_region


def parse_positions(positions_str: str) -> list[int]:
    """Parse comma-separated HxB2 AA position string into sorted ints.

    >>> parse_positions("332")
    [332]
    """
    tokens = [t.strip() for t in positions_str.split(",")]
    positions: list[int] = []
    for token in tokens:
        if not token:
            raise ValueError("Empty position token in positions string")
        try:
            val = int(token)
        except ValueError:
            raise ValueError(f"Non-integer position: {token!r}") from None
        if val <= 0:
            raise ValueError(f"Position must be positive: {val}")
        positions.append(val)
    return sorted(positions)


def parse_regions(region_str: str) -> list[str]:
    """Parse and validate comma-separated Env region names.

    >>> parse_regions(" C1 , V3 ")
    ['C1', 'V3']
    """
    valid_names = {name for name, _, _ in ENV_REGIONS}
    tokens = [t.strip() for t in region_str.split(",")]
    result: list[str] = []
    for token in tokens:
        if not token:
            raise ValueError("Empty region name token")
        if token not in valid_names:
            raise ValueError(
                f"Unknown Env region: {token!r}. "
                f"Valid regions: {sorted(valid_names)}"
            )
        result.append(token)
    return result


def resolve_positions_from_regions(
    regions: list[str],
    reverse_map: dict[int, int],
) -> list[int]:
    """Collect HxB2 positions present in the alignment for given regions.

    >>> rmap = {296: 300, 297: 301, 330: 340, 131: 140, 132: 141}
    >>> resolve_positions_from_regions(["V3"], rmap)
    [296, 297, 330]
    """
    region_bounds = {name: (start, end) for name, start, end in ENV_REGIONS}
    positions: list[int] = []
    for region in regions:
        start, end = region_bounds[region]
        for hxb2_pos in range(start, end + 1):
            if hxb2_pos in reverse_map:
                positions.append(hxb2_pos)
    return sorted(set(positions))


def read_selection_tsv(tsv_path: Path, p_threshold: float) -> list[int]:
    """Read selection_summary.tsv and return significant HxB2 positions.

    Filters to FEL positive sites and MEME episodic sites with
    p_value <= p_threshold.
    """
    positions: set[int] = set()
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            method = row.get("method", "")
            selection_type = row.get("selection_type", "")
            p_value_str = row.get("p_value", "")
            hxb2_str = row.get("hxb2_aa_pos", "")

            if not hxb2_str or hxb2_str.strip() == "":
                continue

            is_fel_positive = method == "FEL" and selection_type == "positive"
            is_meme_episodic = method == "MEME" and selection_type == "episodic"
            if not (is_fel_positive or is_meme_episodic):
                continue

            try:
                p_value = float(p_value_str)
            except (ValueError, TypeError):
                continue
            if p_value > p_threshold:
                continue

            try:
                positions.add(int(hxb2_str))
            except ValueError:
                continue

    return sorted(positions)


def compute_gap_fraction(sequences: list[str], col: int) -> float:
    """Compute fraction of sequences with a gap at the given column.

    >>> compute_gap_fraction([], 0)
    1.0
    """
    if not sequences:
        return 1.0
    n_gaps = sum(1 for seq in sequences if seq[col] in ("-", "."))
    return n_gaps / len(sequences)


def filter_positions_by_gaps(
    hxb2_positions: list[int],
    reverse_map: dict[int, int],
    sample_seqs: list[str],
    max_gap_fraction: float = MAX_GAP_FRACTION,
) -> tuple[list[int], list[int]]:
    """Filter HxB2 positions by gap fraction.

    Returns:
        (selected_cols, valid_hxb2_positions) -- parallel lists of
        alignment column indices and their HxB2 positions.
    """
    selected_cols: list[int] = []
    valid_hxb2: list[int] = []
    for pos in hxb2_positions:
        col = reverse_map[pos]
        gap_frac = compute_gap_fraction(sample_seqs, col)
        if gap_frac > max_gap_fraction:
            continue
        selected_cols.append(col)
        valid_hxb2.append(pos)
    return selected_cols, valid_hxb2


def build_logo_matrix(
    sequences: list[str],
    selected_cols: list[int],
    hxb2_map: list[HxB2Position],
    matrix_type: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Build a logomaker matrix from selected alignment columns.

    Returns:
        (matrix DataFrame, list of x-axis label strings).
    """
    subsequences: list[str] = []
    for seq in sequences:
        subseq = "".join(seq[col] for col in selected_cols)
        subsequences.append(subseq)

    labels: list[str] = []
    for col in selected_cols:
        if col < len(hxb2_map) and hxb2_map[col].hxb2_aa_pos is not None:
            labels.append(str(hxb2_map[col].hxb2_aa_pos))
        else:
            labels.append(f"ins{col}")

    to_type = matrix_type if matrix_type != "counts" else "counts"
    matrix = logomaker.alignment_to_matrix(
        subsequences,
        to_type=to_type,
        characters_to_ignore=".-",
        pseudocount=1.0 if matrix_type == "information" else 0.0,
    )
    matrix.index = range(len(matrix))

    return matrix, labels


def build_region_groups(valid_hxb2: list[int]) -> dict[str, list[int]]:
    """Group matrix column indices by Env region for multi-panel layout."""
    region_groups: dict[str, list[int]] = defaultdict(list)
    for matrix_idx, hxb2_pos in enumerate(valid_hxb2):
        region = get_env_region(hxb2_pos)
        region_name = region if region else "other"
        region_groups[region_name].append(matrix_idx)
    return dict(region_groups)


def find_variant_positions(
    alignment: dict[str, str],
    ref_id: str,
    hxb2_map: list[HxB2Position],
    animal_groups: dict[str, list[str]],
    min_freq: float = 0.50,
) -> list[int]:
    """Find HxB2 positions where any animal has >min_freq variant frequency.

    For each alignment column, checks if any animal's non-gap sequences
    differ from the parental reference above the threshold.

    >>> find_variant_positions({}, "ref", [], {})
    []
    """
    ref_seq = alignment.get(ref_id)
    if ref_seq is None:
        return []

    variant_positions: set[int] = set()
    for col in range(len(ref_seq)):
        ref_aa = ref_seq[col]
        if ref_aa in ("-", "."):
            continue
        if col >= len(hxb2_map) or hxb2_map[col].hxb2_aa_pos is None:
            continue
        hxb2_pos = hxb2_map[col].hxb2_aa_pos
        for _animal, seqs in animal_groups.items():
            non_gap = [s[col] for s in seqs if s[col] not in ("-", ".")]
            if not non_gap:
                continue
            n_diff = sum(1 for aa in non_gap if aa != ref_aa)
            if n_diff / len(non_gap) >= min_freq:
                variant_positions.add(hxb2_pos)
                break

    return sorted(variant_positions)


def resolve_window_cols(
    center_hxb2: int,
    radius: int,
    reverse_map: dict[int, int],
    hxb2_map: list[HxB2Position],
    sample_seqs: list[str],
) -> tuple[list[int], list[str]]:
    """Resolve alignment columns for a window around a center HxB2 position.

    Expands +/- radius alignment columns around the center so insertions
    are captured. Skips columns where gap fraction > MAX_GAP_FRACTION.

    Returns:
        (selected_cols, labels)
    """
    if center_hxb2 not in reverse_map:
        return [], []

    center_col = reverse_map[center_hxb2]
    aln_len = len(sample_seqs[0]) if sample_seqs else 0
    start = max(0, center_col - radius)
    end = min(aln_len, center_col + radius + 1)

    selected_cols: list[int] = []
    labels: list[str] = []
    for col in range(start, end):
        gap_frac = compute_gap_fraction(sample_seqs, col)
        if gap_frac > MAX_GAP_FRACTION:
            continue
        selected_cols.append(col)
        if col < len(hxb2_map) and hxb2_map[col].hxb2_aa_pos is not None:
            labels.append(str(hxb2_map[col].hxb2_aa_pos))
        else:
            labels.append(f"ins{col}")
    return selected_cols, labels
