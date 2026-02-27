"""Click CLI for tlogo."""

import click

from tlogo import __version__


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="tlogo")
def main():
    """Sequence logo generator for aligned FASTA files."""


@main.command()
@click.argument("name", default="World")
def hello(name):
    """Say hello (example command — replace me)."""
    click.echo(f"Hello, {name}!")
