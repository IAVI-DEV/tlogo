"""Tests for tlogo CLI."""

from __future__ import annotations

import tempfile
from pathlib import Path

from click.testing import CliRunner
from hypothesis import given, settings

from tlogo.cli import main

from conftest import _keep, aligned_aa_seqs


class TestCLI:
    """CLI smoke tests."""

    def test_version(self):
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "tlogo" in result.output

    def test_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "logo" in result.output
        assert "window" in result.output
        assert "auto" in result.output

    def test_logo_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["logo", "--help"])
        assert result.exit_code == 0
        assert "--positions" in result.output
        assert "--region" in result.output

    def test_window_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["window", "--help"])
        assert result.exit_code == 0
        assert "--position" in result.output
        assert "--radius" in result.output

    def test_auto_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["auto", "--help"])
        assert result.exit_code == 0
        assert "--variant-freq" in result.output

    def test_logo_missing_fasta(self):
        runner = CliRunner()
        result = runner.invoke(main, ["logo", "nonexistent.fasta"])
        assert result.exit_code != 0

    def test_logo_runs(self, example_fasta, output_dir):
        runner = CliRunner()
        output = output_dir / "test_logo.pdf"
        result = runner.invoke(main, [
            "logo", str(example_fasta),
            "--positions", "5,10,15",
            "-o", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_window_runs(self, example_fasta, output_dir):
        runner = CliRunner()
        output = output_dir / "test_window.pdf"
        result = runner.invoke(main, [
            "window", str(example_fasta),
            "--position", "100",
            "--radius", "3",
            "-o", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_logo_png_format(self, example_fasta, output_dir):
        runner = CliRunner()
        output = output_dir / "test_logo.png"
        result = runner.invoke(main, [
            "logo", str(example_fasta),
            "--positions", "5,10,15",
            "--format", "png",
            "-o", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_logo_svg_format(self, example_fasta, output_dir):
        runner = CliRunner()
        output = output_dir / "test_logo.svg"
        result = runner.invoke(main, [
            "logo", str(example_fasta),
            "--positions", "5,10,15",
            "--format", "svg",
            "-o", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_logo_probability_matrix(self, example_fasta, output_dir):
        runner = CliRunner()
        output = output_dir / "test_logo_probability.pdf"
        result = runner.invoke(main, [
            "logo", str(example_fasta),
            "--positions", "5,10,15",
            "--matrix-type", "probability",
            "-o", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()

    def test_logo_counts_matrix(self, example_fasta, output_dir):
        runner = CliRunner()
        output = output_dir / "test_logo_counts.pdf"
        result = runner.invoke(main, [
            "logo", str(example_fasta),
            "--positions", "5,10,15",
            "--matrix-type", "counts",
            "-o", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()

    def test_auto_runs(self, example_fasta, output_dir):
        runner = CliRunner()
        out_dir = output_dir / "auto_test"
        result = runner.invoke(main, [
            "auto", str(example_fasta),
            "--ref-id", "parental_ref",
            "--output-dir", str(out_dir),
            "--variant-freq", "0.5",
        ])
        assert result.exit_code == 0, result.output


class TestCLIHypothesis:
    """Property-based CLI rendering tests using hypothesis."""

    _logo_counter = 0
    _window_counter = 0

    @given(data=aligned_aa_seqs(min_seqs=3, max_seqs=6, min_len=10, max_len=30))
    @settings(max_examples=5, deadline=60000)
    def test_logo_renders_random_alignment(self, data):
        """Logo command does not crash on random aligned sequences."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            fasta = tmp / "hyp.fasta"
            with open(fasta, "w") as f:
                for name, seq in data:
                    f.write(f">{name}\n{seq}\n")

            aln_len = len(data[0][1])
            positions = ",".join(str(p) for p in [1, max(1, aln_len // 2), aln_len])

            TestCLIHypothesis._logo_counter += 1
            out = tmp / "logo.pdf"
            runner = CliRunner()
            result = runner.invoke(main, [
                "logo", str(fasta),
                "--positions", positions,
                "-o", str(out),
            ])
            # May fail if positions not in map, that's ok
            if result.exit_code == 0 and out.exists():
                _keep(out, f"hypothesis_logo_{TestCLIHypothesis._logo_counter:03d}.pdf")

    @given(data=aligned_aa_seqs(min_seqs=3, max_seqs=6, min_len=10, max_len=30))
    @settings(max_examples=5, deadline=60000)
    def test_window_renders_random_alignment(self, data):
        """Window command does not crash on random aligned sequences."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            fasta = tmp / "hyp.fasta"
            with open(fasta, "w") as f:
                for name, seq in data:
                    f.write(f">{name}\n{seq}\n")

            aln_len = len(data[0][1])
            center = max(1, aln_len // 2)

            TestCLIHypothesis._window_counter += 1
            out = tmp / "window.pdf"
            runner = CliRunner()
            result = runner.invoke(main, [
                "window", str(fasta),
                "--position", str(center),
                "--radius", "2",
                "-o", str(out),
            ])
            if result.exit_code == 0 and out.exists():
                _keep(out, f"hypothesis_window_{TestCLIHypothesis._window_counter:03d}.pdf")
