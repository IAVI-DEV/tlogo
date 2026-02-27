"""Tests for tlogo.hxb2 module."""

from __future__ import annotations

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from tlogo.hxb2 import (
    HxB2Position,
    build_hxb2_map,
    build_reverse_hxb2_map,
    get_env_region,
)

from conftest import AA_ALPHABET, aligned_aa_seqs


class TestGetEnvRegion:
    def test_signal_peptide(self):
        assert get_env_region(1) == "SP"
        assert get_env_region(30) == "SP"

    def test_v_loops(self):
        assert get_env_region(131) == "V1"
        assert get_env_region(296) == "V3"
        assert get_env_region(460) == "V5"

    def test_constant_regions(self):
        assert get_env_region(31) == "C1"
        assert get_env_region(197) == "C2"

    def test_gp41(self):
        assert get_env_region(512) == "gp41"
        assert get_env_region(856) == "gp41"

    def test_out_of_range(self):
        assert get_env_region(0) is None
        assert get_env_region(857) is None
        assert get_env_region(-1) is None

    @given(pos=st.integers(min_value=1, max_value=856))
    @settings(max_examples=50)
    def test_all_valid_positions_return_region(self, pos):
        """Every position 1-856 maps to some region."""
        region = get_env_region(pos)
        assert region is not None

    @given(pos=st.integers(min_value=857, max_value=10000))
    @settings(max_examples=20)
    def test_beyond_gp41_returns_none(self, pos):
        assert get_env_region(pos) is None


class TestBuildHxB2Map:
    def test_simple_alignment(self):
        alignment = {
            "HxB2": "ABC",
            "sample": "ABD",
        }
        hmap = build_hxb2_map(alignment)
        assert len(hmap) == 3
        assert hmap[0].hxb2_aa_pos == 1
        assert hmap[1].hxb2_aa_pos == 2
        assert hmap[2].hxb2_aa_pos == 3

    def test_gap_in_hxb2(self):
        alignment = {
            "HxB2": "A-B",
            "sample": "ACB",
        }
        hmap = build_hxb2_map(alignment)
        assert hmap[0].hxb2_aa_pos == 1
        assert hmap[1].hxb2_aa_pos is None  # insertion column
        assert hmap[1].hxb2_residue == "-"
        assert hmap[2].hxb2_aa_pos == 2

    def test_missing_hxb2_raises(self):
        alignment = {"sample": "ABC"}
        with pytest.raises(ValueError, match="not found"):
            build_hxb2_map(alignment)

    def test_custom_hxb2_id(self):
        alignment = {"MyRef": "ABC", "sample": "ABD"}
        hmap = build_hxb2_map(alignment, hxb2_id="MyRef")
        assert len(hmap) == 3

    @given(data=aligned_aa_seqs(min_seqs=2, max_seqs=5, min_len=3, max_len=30))
    @settings(max_examples=30)
    def test_map_length_matches_alignment(self, data):
        """HxB2 map always has one entry per alignment column."""
        alignment = {name: seq for name, seq in data}
        hmap = build_hxb2_map(alignment)
        aln_len = len(data[0][1])
        assert len(hmap) == aln_len

    @given(data=aligned_aa_seqs(min_seqs=2, max_seqs=5, min_len=3, max_len=30))
    @settings(max_examples=30)
    def test_positions_are_sequential(self, data):
        """Non-None hxb2_aa_pos values are strictly sequential 1..N."""
        alignment = {name: seq for name, seq in data}
        hmap = build_hxb2_map(alignment)
        positions = [p.hxb2_aa_pos for p in hmap if p.hxb2_aa_pos is not None]
        assert positions == list(range(1, len(positions) + 1))

    @given(data=aligned_aa_seqs(min_seqs=2, max_seqs=5, min_len=3, max_len=30))
    @settings(max_examples=30)
    def test_non_gap_count_matches_hxb2_length(self, data):
        """Number of non-None positions equals non-gap residues in HxB2."""
        alignment = {name: seq for name, seq in data}
        hxb2_seq = alignment["HxB2"]
        hmap = build_hxb2_map(alignment)
        non_gap = sum(1 for c in hxb2_seq if c != "-")
        mapped = sum(1 for p in hmap if p.hxb2_aa_pos is not None)
        assert mapped == non_gap


class TestBuildReverseHxB2Map:
    def test_basic(self):
        positions = [
            HxB2Position(0, 10, "C1", "A"),
            HxB2Position(1, None, None, "-"),
            HxB2Position(2, 11, "C1", "B"),
        ]
        rmap = build_reverse_hxb2_map(positions)
        assert rmap == {10: 0, 11: 2}

    def test_empty(self):
        assert build_reverse_hxb2_map([]) == {}

    @given(data=aligned_aa_seqs(min_seqs=2, max_seqs=4, min_len=3, max_len=30))
    @settings(max_examples=30)
    def test_roundtrip_map_reverse(self, data):
        """Reverse map correctly inverts forward map for non-None entries."""
        alignment = {name: seq for name, seq in data}
        hmap = build_hxb2_map(alignment)
        rmap = build_reverse_hxb2_map(hmap)
        for p in hmap:
            if p.hxb2_aa_pos is not None:
                assert rmap[p.hxb2_aa_pos] == p.alignment_col
