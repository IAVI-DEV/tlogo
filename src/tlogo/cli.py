"""Click CLI for tlogo."""

from __future__ import annotations

from pathlib import Path

import click

from tlogo import __version__


# ---------------------------------------------------------------------------
# Shared options
# ---------------------------------------------------------------------------

def _common_options(func):
    """Stack shared Click options onto a command."""
    func = click.option(
        "--color-scheme",
        type=click.Choice(["chemistry", "charge", "hydrophobicity"]),
        default="chemistry", show_default=True,
        help="Amino acid color scheme.",
    )(func)
    func = click.option(
        "--matrix-type",
        type=click.Choice(["information", "probability", "counts"]),
        default="information", show_default=True,
        help="Type of matrix for logo height.",
    )(func)
    func = click.option(
        "--ref-id", default=None,
        help="Parental reference seq ID (auto-detected if omitted).",
    )(func)
    func = click.option(
        "--hxb2-id", default="HxB2", show_default=True,
        help="HxB2 seq ID in the alignment.",
    )(func)
    func = click.option(
        "--delimiter", default="_", show_default=True,
        help="Seq ID field separator for animal name extraction.",
    )(func)
    func = click.option(
        "--field", type=int, default=0, show_default=True,
        help="0-based field index for animal name.",
    )(func)
    func = click.option(
        "--format", "fmt",
        type=click.Choice(["pdf", "png", "svg"]),
        default="pdf", show_default=True,
        help="Output format.",
    )(func)
    return func


# ---------------------------------------------------------------------------
# Group
# ---------------------------------------------------------------------------

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="tlogo")
def main() -> None:
    """tlogo -- Sequence logo plots for aligned FASTA files (HIV Env / HxB2)."""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_alignment(fasta: Path, hxb2_id: str, ref_id: str | None):
    """Parse alignment, build HxB2 map, detect ref, collect sample seqs."""
    from tlogo.hxb2 import build_hxb2_map
    from tlogo.io import detect_reference_id, group_sequences_by_animal, parse_alignment

    alignment = parse_alignment(fasta)
    if hxb2_id not in alignment:
        raise click.BadParameter(
            f"HxB2 sequence '{hxb2_id}' not found in {fasta}",
            param_hint="--hxb2-id",
        )

    hxb2_map = build_hxb2_map(alignment, hxb2_id)
    detected_ref = detect_reference_id(alignment, ref_id, hxb2_id)

    return alignment, hxb2_map, detected_ref


def _collect_samples(
    alignment: dict[str, str],
    hxb2_id: str,
    ref_id: str | None,
    delimiter: str,
    field: int,
):
    """Build exclude set, sample list, and animal groups."""
    from tlogo.io import group_sequences_by_animal

    exclude = {hxb2_id}
    if ref_id is not None:
        exclude.add(ref_id)

    animal_groups = group_sequences_by_animal(
        alignment, exclude, delimiter, field,
    )
    sample_seqs = [s for seqs in animal_groups.values() for s in seqs]

    return exclude, sample_seqs, animal_groups


# ---------------------------------------------------------------------------
# tlogo logo
# ---------------------------------------------------------------------------

@main.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--positions", default=None,
    help="Comma-separated HxB2 AA positions (e.g. 130,160,332).",
)
@click.option(
    "--region", default=None,
    help="Env region(s) (e.g. V3 or V1,V2,V3).",
)
@click.option(
    "--selection-tsv", type=click.Path(exists=True, path_type=Path),
    default=None,
    help="selection_summary.tsv for automatic position selection.",
)
@click.option(
    "--p-threshold", type=float, default=0.1, show_default=True,
    help="P-value threshold for --selection-tsv.",
)
@click.option(
    "--max-gap-fraction", type=float, default=0.9, show_default=True,
    help="Skip positions above this gap fraction.",
)
@click.option("--output", "-o", type=click.Path(path_type=Path), default=None)
@_common_options
def logo(
    fasta: Path,
    positions: str | None,
    region: str | None,
    selection_tsv: Path | None,
    p_threshold: float,
    max_gap_fraction: float,
    output: Path | None,
    color_scheme: str,
    matrix_type: str,
    ref_id: str | None,
    hxb2_id: str,
    delimiter: str,
    field: int,
    fmt: str,
) -> None:
    """Standard logo plot at selected HxB2 positions.

    FASTA is a protein alignment containing HxB2 for coordinate mapping.
    """
    from tlogo.logo import process_standard_logo
    from tlogo.matrix import parse_positions, parse_regions, read_selection_tsv

    alignment, hxb2_map, detected_ref = _load_alignment(fasta, hxb2_id, ref_id)
    exclude, sample_seqs, animal_groups = _collect_samples(
        alignment, hxb2_id, detected_ref, delimiter, field,
    )

    if not sample_seqs:
        click.echo("No sample sequences found.", err=True)
        raise SystemExit(1)

    # Resolve positions
    resolved_positions: list[int] | None = None
    resolved_regions: list[str] | None = None

    if positions is not None:
        resolved_positions = parse_positions(positions)
    elif selection_tsv is not None:
        resolved_positions = read_selection_tsv(selection_tsv, p_threshold)
    elif region is not None:
        resolved_regions = parse_regions(region)

    # Default output
    if output is None:
        output = fasta.parent / f"logo_plot.{fmt}"

    title = fasta.stem
    ok = process_standard_logo(
        alignment, hxb2_map, sample_seqs, output,
        title=title,
        positions=resolved_positions,
        regions=resolved_regions,
        matrix_type=matrix_type,
        color_scheme=color_scheme,
        max_gap_fraction=max_gap_fraction,
    )

    if ok:
        click.echo(f"Output: {output}")
    else:
        click.echo("No positions passed filters.", err=True)
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# tlogo window
# ---------------------------------------------------------------------------

@main.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--position", "-p", type=int, required=True,
    help="Center HxB2 AA position for the window.",
)
@click.option(
    "--radius", "-r", type=int, default=4, show_default=True,
    help="Window radius in alignment columns.",
)
@click.option("--output", "-o", type=click.Path(path_type=Path), default=None)
@_common_options
def window(
    fasta: Path,
    position: int,
    radius: int,
    output: Path | None,
    color_scheme: str,
    matrix_type: str,
    ref_id: str | None,
    hxb2_id: str,
    delimiter: str,
    field: int,
    fmt: str,
) -> None:
    """Per-animal window logo at a focused HxB2 position +/- radius.

    FASTA is a protein alignment containing HxB2 for coordinate mapping.
    """
    from tlogo.window import process_window_mode

    alignment, hxb2_map, detected_ref = _load_alignment(fasta, hxb2_id, ref_id)
    exclude, sample_seqs, animal_groups = _collect_samples(
        alignment, hxb2_id, detected_ref, delimiter, field,
    )

    if not sample_seqs:
        click.echo("No sample sequences found.", err=True)
        raise SystemExit(1)

    if output is None:
        output = fasta.parent / f"logo_window_{position}.{fmt}"

    title = fasta.stem
    ok = process_window_mode(
        alignment, hxb2_map, position, radius, animal_groups,
        output, title, matrix_type, color_scheme,
    )

    if ok:
        click.echo(f"Output: {output}")
    else:
        click.echo(
            f"HxB2 position {position} not found or all columns filtered.",
            err=True,
        )
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# tlogo auto
# ---------------------------------------------------------------------------

@main.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--radius", "-r", type=int, default=4, show_default=True,
    help="Window radius in alignment columns.",
)
@click.option(
    "--variant-freq", type=float, default=0.50, show_default=True,
    help="Minimum per-animal variant frequency threshold.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(path_type=Path), default=None,
    help="Output directory for window logo files.",
)
@_common_options
def auto(
    fasta: Path,
    radius: int,
    variant_freq: float,
    output_dir: Path | None,
    color_scheme: str,
    matrix_type: str,
    ref_id: str | None,
    hxb2_id: str,
    delimiter: str,
    field: int,
    fmt: str,
) -> None:
    """Auto-detect variant positions and generate window logos for each.

    Scans all alignment positions for per-animal variant frequency above
    the threshold, then generates a window logo plot for each.

    FASTA is a protein alignment containing HxB2 for coordinate mapping.
    """
    from tlogo.window import process_auto_variants

    alignment, hxb2_map, detected_ref = _load_alignment(fasta, hxb2_id, ref_id)

    if detected_ref is None:
        click.echo(
            "Could not detect parental reference. "
            "Use --ref-id to specify it.",
            err=True,
        )
        raise SystemExit(1)

    exclude, sample_seqs, animal_groups = _collect_samples(
        alignment, hxb2_id, detected_ref, delimiter, field,
    )

    if not sample_seqs:
        click.echo("No sample sequences found.", err=True)
        raise SystemExit(1)

    if output_dir is None:
        output_dir = fasta.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    title = fasta.stem
    generated = process_auto_variants(
        alignment, hxb2_map, animal_groups, detected_ref,
        output_dir, title, radius, variant_freq,
        matrix_type, color_scheme, fmt=fmt,
    )

    if generated:
        click.echo(
            f"Generated {len(generated)} window logos at HxB2 positions: "
            f"{generated}"
        )
        click.echo(f"Output directory: {output_dir}")
    else:
        click.echo("No variant positions detected.", err=True)
