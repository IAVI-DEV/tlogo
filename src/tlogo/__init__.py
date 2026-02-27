"""tlogo — Sequence logo generator for aligned FASTA files."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tlogo")
except PackageNotFoundError:
    __version__ = "0.0.0"
