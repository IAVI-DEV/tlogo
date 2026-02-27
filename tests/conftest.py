"""Shared fixtures for tlogo tests."""

from __future__ import annotations

from pathlib import Path

import pytest

OUTPUT_DIR = Path(__file__).parent / "output"


@pytest.fixture
def output_dir() -> Path:
    OUTPUT_DIR.mkdir(exist_ok=True)
    return OUTPUT_DIR


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
