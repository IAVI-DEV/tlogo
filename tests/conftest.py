"""Shared fixtures and hypothesis strategies for tlogo tests."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from hypothesis import strategies as st

OUTPUT_DIR = Path(__file__).parent / "output"
EXAMPLES_DIR = Path(__file__).parent.parent / "examples"

# Amino acid alphabet (standard 20 + gap)
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_WITH_GAP = AA_ALPHABET + "-"


# ---------------------------------------------------------------------------
# Output directory fixtures (tview pattern)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session", autouse=True)
def ensure_output_dir():
    """Create tests/output/ once per session."""
    OUTPUT_DIR.mkdir(exist_ok=True)


@pytest.fixture
def output_dir() -> Path:
    """Persistent output directory for visual inspection."""
    return OUTPUT_DIR


def _keep(src: str | Path, name: str) -> None:
    """Copy rendered file into tests/output/ for visual inspection."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    shutil.copy2(str(src), OUTPUT_DIR / name)


# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

seq_name = st.text(
    alphabet="abcdefghijklmnopqrstuvwxyz0123456789",
    min_size=1,
    max_size=12,
)


@st.composite
def aligned_aa_seqs(
    draw,
    min_seqs: int = 2,
    max_seqs: int = 8,
    min_len: int = 5,
    max_len: int = 50,
    gap_freq: float = 0.1,
):
    """Generate aligned AA sequences (same length, with optional gaps).

    Returns list of (name, sequence) tuples with an HxB2 entry first.
    """
    length = draw(st.integers(min_value=min_len, max_value=max_len))
    n_seqs = draw(st.integers(min_value=min_seqs, max_value=max_seqs))

    # Build HxB2 reference (no gaps in HxB2 for valid mapping)
    hxb2_seq = draw(
        st.text(alphabet=AA_ALPHABET, min_size=length, max_size=length)
    )
    seqs = [("HxB2", hxb2_seq)]

    # Build sample sequences with animal-style names
    for i in range(n_seqs - 1):
        animal = draw(seq_name)
        seq = draw(
            st.text(alphabet=AA_WITH_GAP, min_size=length, max_size=length)
        )
        seqs.append((f"{animal}_wk{i}_00{i}", seq))

    return seqs


@st.composite
def hiv_like_alignment(
    draw,
    min_animals: int = 1,
    max_animals: int = 4,
    seqs_per_animal: int = 2,
    min_len: int = 10,
    max_len: int = 30,
):
    """Generate an alignment dict with HxB2, optional ref, and animal groups.

    Returns (alignment_dict, has_ref) tuple.
    """
    length = draw(st.integers(min_value=min_len, max_value=max_len))
    n_animals = draw(st.integers(min_value=min_animals, max_value=max_animals))
    has_ref = draw(st.booleans())

    alignment: dict[str, str] = {}

    # HxB2 reference (no gaps)
    hxb2_seq = draw(
        st.text(alphabet=AA_ALPHABET, min_size=length, max_size=length)
    )
    alignment["HxB2"] = hxb2_seq

    # Optional parental reference
    if has_ref:
        ref_seq = draw(
            st.text(alphabet=AA_WITH_GAP, min_size=length, max_size=length)
        )
        alignment["parental_ref"] = ref_seq

    # Animal sequences
    for a in range(n_animals):
        animal_prefix = f"rh{1000 + a}"
        for s in range(seqs_per_animal):
            seq = draw(
                st.text(
                    alphabet=AA_WITH_GAP, min_size=length, max_size=length
                )
            )
            alignment[f"{animal_prefix}_wk4_{s:03d}"] = seq

    return alignment, has_ref


# ---------------------------------------------------------------------------
# File-based fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def write_fasta(tmp_path):
    """Fixture that returns a helper to write FASTA files."""

    def _write(seqs: list[tuple[str, str]], name: str = "test.fasta") -> Path:
        path = tmp_path / name
        with open(path, "w", encoding="utf-8") as fh:
            for seq_name, seq in seqs:
                fh.write(f">{seq_name}\n{seq}\n")
        return path

    return _write


# A minimal HIV-like protein alignment for testing
_MINI_SEQS = [
    ("HxB2",              "MRVKETKLWVTVYYGVPVW"),
    ("parental_ref",      "MRVKETKLWVTVYYGVPVW"),
    ("animalA_wk8_001",   "MRVKETKLWVTVYYGVPVW"),
    ("animalA_wk8_002",   "MRVRETKLWVTVYYGVPVW"),
    ("animalB_wk4_001",   "MRVKETRLWVTVYYGVPVW"),
    ("animalB_wk4_002",   "MRVKETRLWVTVYYGVPVW"),
]


@pytest.fixture
def mini_alignment(write_fasta) -> Path:
    """A minimal HIV-like protein alignment with HxB2 and animal sequences."""
    return write_fasta(_MINI_SEQS, "mini_alignment.fasta")


@pytest.fixture
def mini_alignment_dict() -> dict[str, str]:
    """The mini alignment as a pre-parsed dict."""
    return {name: seq for name, seq in _MINI_SEQS}


@pytest.fixture
def example_fasta() -> Path:
    """Path to the example.fasta in the examples directory."""
    return EXAMPLES_DIR / "example.fasta"
