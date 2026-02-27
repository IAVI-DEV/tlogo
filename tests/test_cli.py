"""Tests for tlogo CLI."""

from click.testing import CliRunner

from tlogo.cli import main


class TestCLI:
    def test_version(self):
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "tlogo" in result.output

    def test_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "Sequence logo generator for aligned FASTA files" in result.output

    def test_hello_default(self):
        runner = CliRunner()
        result = runner.invoke(main, ["hello"])
        assert result.exit_code == 0
        assert "Hello, World!" in result.output

    def test_hello_with_name(self):
        runner = CliRunner()
        result = runner.invoke(main, ["hello", "Bioinformatics"])
        assert result.exit_code == 0
        assert "Hello, Bioinformatics!" in result.output
