"""Tests for tlogo.matrix module."""

from __future__ import annotations

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from tlogo.hxb2 import HxB2Position, build_hxb2_map, build_reverse_hxb2_map
from tlogo.matrix import (
    build_region_groups,
    compute_gap_fraction,
    filter_positions_by_gaps,
    find_variant_positions,
    parse_positions,
    parse_regions,
    resolve_positions_from_regions,
    resolve_window_cols,
)

from conftest import AA_ALPHABET, AA_WITH_GAP, aligned_aa_seqs


class TestParsePositions:
    def test_basic(self):
        assert parse_positions("130,160,332") == [130, 160, 332]

    def test_whitespace(self):
        assert parse_positions(" 50 , 100 , 200 ") == [50, 100, 200]

    def test_single(self):
        assert parse_positions("332") == [332]

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="positive"):
            parse_positions("-5")

    def test_non_integer_raises(self):
        with pytest.raises(ValueError, match="Non-integer"):
            parse_positions("abc")

    @given(
        positions=st.lists(
            st.integers(min_value=1, max_value=856),
            min_size=1,
            max_size=10,
        )
    )
    @settings(max_examples=30)
    def test_roundtrip_parse(self, positions):
        """Parsing comma-joined positive ints always returns sorted list."""
        s = ",".join(str(p) for p in positions)
        result = parse_positions(s)
        assert result == sorted(positions)

    @given(
        positions=st.lists(
            st.integers(min_value=1, max_value=856),
            min_size=1,
            max_size=10,
        )
    )
    @settings(max_examples=20)
    def test_output_is_sorted(self, positions):
        """Output is always sorted regardless of input order."""
        s = ",".join(str(p) for p in positions)
        result = parse_positions(s)
        assert result == sorted(result)


class TestParseRegions:
    def test_single(self):
        assert parse_regions("V3") == ["V3"]

    def test_multiple(self):
        assert parse_regions("V1,V2,V3") == ["V1", "V2", "V3"]

    def test_whitespace(self):
        assert parse_regions(" C1 , V3 ") == ["C1", "V3"]

    def test_invalid_raises(self):
        with pytest.raises(ValueError, match="Unknown"):
            parse_regions("X99")


class TestResolvePositionsFromRegions:
    def test_v3(self):
        rmap = {296: 300, 297: 301, 330: 340, 131: 140}
        result = resolve_positions_from_regions(["V3"], rmap)
        assert result == [296, 297, 330]

    def test_multiple_regions(self):
        rmap = {131: 0, 296: 1, 297: 2}
        result = resolve_positions_from_regions(["V1", "V3"], rmap)
        assert 131 in result
        assert 296 in result


class TestComputeGapFraction:
    def test_all_gaps(self):
        assert compute_gap_fraction(["---", "---", "---"], 0) == 1.0

    def test_no_gaps(self):
        assert compute_gap_fraction(["ABC", "DEF", "GHI"], 0) == 0.0

    def test_mixed(self):
        assert compute_gap_fraction(["A-C", "-BC", "ABC", "A-C"], 1) == 0.5

    def test_empty(self):
        assert compute_gap_fraction([], 0) == 1.0

    @given(
        n_seqs=st.integers(min_value=1, max_value=20),
        data=st.data(),
    )
    @settings(max_examples=30)
    def test_fraction_in_range(self, n_seqs, data):
        """Gap fraction is always between 0.0 and 1.0."""
        seqs = [
            data.draw(st.text(alphabet=AA_WITH_GAP, min_size=1, max_size=1))
            for _ in range(n_seqs)
        ]
        frac = compute_gap_fraction(seqs, 0)
        assert 0.0 <= frac <= 1.0

    @given(n=st.integers(min_value=1, max_value=20))
    @settings(max_examples=20)
    def test_all_aa_gives_zero(self, n):
        """All-AA column has gap fraction 0."""
        seqs = ["A" for _ in range(n)]
        assert compute_gap_fraction(seqs, 0) == 0.0

    @given(n=st.integers(min_value=1, max_value=20))
    @settings(max_examples=20)
    def test_all_gaps_gives_one(self, n):
        """All-gap column has gap fraction 1."""
        seqs = ["-" for _ in range(n)]
        assert compute_gap_fraction(seqs, 0) == 1.0


class TestFilterPositionsByGaps:
    def test_filters_high_gaps(self):
        rmap = {10: 0, 11: 1}
        seqs = ["A-", "A-", "A-"]  # col 1 is 100% gaps
        cols, hxb2 = filter_positions_by_gaps([10, 11], rmap, seqs)
        assert 10 in hxb2
        assert 11 not in hxb2


class TestBuildRegionGroups:
    def test_groups_by_region(self):
        groups = build_region_groups([131, 132, 296, 297])
        assert "V1" in groups
        assert "V3" in groups
        assert len(groups["V1"]) == 2
        assert len(groups["V3"]) == 2


class TestResolveWindowCols:
    def test_basic(self):
        hmap = [HxB2Position(i, i + 5, "C1", "A") for i in range(15)]
        rmap = {i + 5: i for i in range(15)}
        seqs = ["A" * 15]
        cols, labels = resolve_window_cols(10, 2, rmap, hmap, seqs)
        assert len(cols) == 5  # 5-2=3 to 5+2=7
        assert "10" in labels


class TestFindVariantPositions:
    def test_no_variants(self):
        alignment = {
            "ref": "ABC",
            "s1": "ABC",
            "s2": "ABC",
        }
        hmap = [
            HxB2Position(0, 1, "SP", "A"),
            HxB2Position(1, 2, "SP", "B"),
            HxB2Position(2, 3, "SP", "C"),
        ]
        groups = {"animal1": ["ABC", "ABC"]}
        result = find_variant_positions(alignment, "ref", hmap, groups)
        assert result == []

    def test_finds_variant(self):
        alignment = {
            "ref": "ABC",
            "s1": "AXC",
            "s2": "AXC",
        }
        hmap = [
            HxB2Position(0, 1, "SP", "A"),
            HxB2Position(1, 2, "SP", "B"),
            HxB2Position(2, 3, "SP", "C"),
        ]
        groups = {"animal1": ["AXC", "AXC"]}
        result = find_variant_positions(alignment, "ref", hmap, groups)
        assert 2 in result  # position 2 has 100% variant

    def test_empty(self):
        assert find_variant_positions({}, "ref", [], {}) == []

    @given(data=aligned_aa_seqs(min_seqs=3, max_seqs=6, min_len=5, max_len=20))
    @settings(max_examples=20)
    def test_variant_positions_are_valid_hxb2(self, data):
        """All returned variant positions are valid HxB2 positions."""
        alignment = {name: seq for name, seq in data}
        hmap = build_hxb2_map(alignment)
        # Use first non-HxB2 seq as ref
        ref_id = data[1][0]
        groups = {"animal1": [seq for name, seq in data[2:]]}
        result = find_variant_positions(alignment, ref_id, hmap, groups)
        valid_positions = {p.hxb2_aa_pos for p in hmap if p.hxb2_aa_pos is not None}
        for pos in result:
            assert pos in valid_positions
