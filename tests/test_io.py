"""Tests for tlogo.io module."""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from tlogo.io import (
    detect_reference_id,
    extract_animal_name,
    group_sequences_by_animal,
    parse_alignment,
    sort_animal_groups,
)

from conftest import hiv_like_alignment


class TestParseAlignment:
    def test_reads_fasta(self, mini_alignment):
        aln = parse_alignment(mini_alignment)
        assert "HxB2" in aln
        assert len(aln) == 6
        # All sequences should be the same length
        lengths = {len(seq) for seq in aln.values()}
        assert len(lengths) == 1


class TestExtractAnimalName:
    def test_default(self):
        assert extract_animal_name("animalA_wk8_001") == "animalA"

    def test_custom_field(self):
        assert extract_animal_name("animalA_wk8_001", field=1) == "wk8"

    def test_no_delimiter(self):
        assert extract_animal_name("singlename") == "singlename"

    def test_custom_delimiter(self):
        assert extract_animal_name("animalA-wk8", delimiter="-") == "animalA"

    @given(
        parts=st.lists(
            st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", min_size=1, max_size=8),
            min_size=1,
            max_size=5,
        )
    )
    @settings(max_examples=30)
    def test_field0_is_first_part(self, parts):
        """Field 0 always returns the first part of a delimited ID."""
        seq_id = "_".join(parts)
        assert extract_animal_name(seq_id) == parts[0]

    @given(
        parts=st.lists(
            st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", min_size=1, max_size=8),
            min_size=2,
            max_size=5,
        ),
        field=st.integers(min_value=0, max_value=4),
    )
    @settings(max_examples=30)
    def test_valid_field_returns_part(self, parts, field):
        """Valid field indices return the correct part."""
        seq_id = "_".join(parts)
        result = extract_animal_name(seq_id, field=field)
        if field < len(parts):
            assert result == parts[field]
        else:
            assert result == seq_id


class TestGroupSequencesByAnimal:
    def test_groups_correctly(self, mini_alignment_dict):
        exclude = {"HxB2", "parental_ref"}
        groups = group_sequences_by_animal(mini_alignment_dict, exclude)
        assert "animalA" in groups
        assert "animalB" in groups
        assert len(groups["animalA"]) == 2
        assert len(groups["animalB"]) == 2

    def test_excludes_ids(self, mini_alignment_dict):
        exclude = {"HxB2", "parental_ref"}
        groups = group_sequences_by_animal(mini_alignment_dict, exclude)
        all_seqs = [s for seqs in groups.values() for s in seqs]
        assert len(all_seqs) == 4  # 6 total - 2 excluded

    @given(data=hiv_like_alignment(min_animals=1, max_animals=3, seqs_per_animal=2))
    @settings(max_examples=20)
    def test_exclude_removes_all_excluded(self, data):
        """Excluded IDs never appear in grouped output."""
        alignment, has_ref = data
        exclude = {"HxB2"}
        if has_ref:
            exclude.add("parental_ref")
        groups = group_sequences_by_animal(alignment, exclude)
        all_ids_in_groups = sum(len(seqs) for seqs in groups.values())
        expected = len(alignment) - len(exclude)
        assert all_ids_in_groups == expected


class TestSortAnimalGroups:
    def test_self_first(self):
        result = sort_animal_groups(
            ["animalA", "parental", "Rec5050", "animalB"], "parental"
        )
        assert result[0] == "parental"

    def test_recombinants_second(self):
        result = sort_animal_groups(["beta", "alpha", "RecX", "RecA"])
        assert result[:2] == ["RecA", "RecX"]

    def test_alphabetical_others(self):
        result = sort_animal_groups(["beta", "alpha", "gamma"])
        assert result == ["alpha", "beta", "gamma"]

    def test_no_self(self):
        result = sort_animal_groups(["beta", "alpha"])
        assert result == ["alpha", "beta"]

    @given(
        names=st.lists(
            st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", min_size=1, max_size=8),
            min_size=1,
            max_size=10,
            unique=True,
        )
    )
    @settings(max_examples=30)
    def test_preserves_all_elements(self, names):
        """Sorting preserves all elements (no loss or duplication)."""
        result = sort_animal_groups(names)
        assert sorted(result) == sorted(names)

    @given(
        names=st.lists(
            st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", min_size=1, max_size=8),
            min_size=2,
            max_size=8,
            unique=True,
        )
    )
    @settings(max_examples=20)
    def test_self_always_first(self, names):
        """When self_name is provided and present, it's always first."""
        self_name = names[0]
        result = sort_animal_groups(names, self_name)
        assert result[0] == self_name


class TestDetectReferenceId:
    def test_explicit(self):
        aln = {"HxB2": "A", "parental_ref": "A", "sample": "A"}
        assert detect_reference_id(aln, "parental_ref") == "parental_ref"

    def test_explicit_missing(self):
        aln = {"HxB2": "A", "sample": "A"}
        assert detect_reference_id(aln, "parental_ref") is None

    def test_auto_detect_ref_suffix(self):
        aln = {"HxB2": "A", "parental_ref": "A", "sample": "A"}
        assert detect_reference_id(aln) == "parental_ref"

    def test_no_ref(self):
        aln = {"HxB2": "A", "sample": "A"}
        assert detect_reference_id(aln) is None

    @given(data=hiv_like_alignment())
    @settings(max_examples=20)
    def test_auto_detect_finds_ref_suffix(self, data):
        """Auto-detection finds _ref suffix when present."""
        alignment, has_ref = data
        result = detect_reference_id(alignment)
        if has_ref:
            assert result == "parental_ref"
        else:
            assert result is None
