"""Standard logo plot rendering (single panel and multi-panel by region)."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import logomaker
import pandas as pd

from tlogo.constants import (
    ENV_REGIONS,
    FONT_AXIS_LABEL,
    FONT_TICK_LABEL,
    FONT_TITLE,
    NATURE_DOUBLE_COL,
    NATURE_SINGLE_COL,
    PANEL_HEIGHT,
    PANEL_HEIGHT_MULTI,
    V_LOOP_REGIONS,
)
from tlogo.hxb2 import HxB2Position, build_hxb2_map, build_reverse_hxb2_map
from tlogo.matrix import (
    build_logo_matrix,
    build_region_groups,
    filter_positions_by_gaps,
    resolve_positions_from_regions,
)


def apply_nature_style() -> None:
    """Apply Nature-journal-compliant matplotlib rcParams."""
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family":         "sans-serif",
        "font.sans-serif":     ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size":           7.0,
        "axes.titlesize":      8.0,
        "axes.labelsize":      7.0,
        "xtick.labelsize":     6.0,
        "ytick.labelsize":     7.0,
        "xtick.direction":     "out",
        "ytick.direction":     "out",
        "xtick.major.size":    3.0,
        "ytick.major.size":    3.0,
        "xtick.major.width":   0.6,
        "ytick.major.width":   0.6,
        "xtick.major.pad":     2.0,
        "ytick.major.pad":     2.0,
        "axes.spines.top":     False,
        "axes.spines.right":   False,
        "axes.linewidth":      0.6,
        "figure.dpi":          150,
        "savefig.dpi":         300,
        "pdf.fonttype":        42,
        "ps.fonttype":         42,
    })


def _y_axis_label(matrix_type: str) -> str:
    """Return y-axis label for the given matrix type."""
    return {
        "information": "Information (bits)",
        "probability": "Frequency",
        "counts": "Count",
    }.get(matrix_type, "Value")


def _save_figure(fig: plt.Figure, output_path: Path, fmt: str) -> None:
    """Save figure to the specified format."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "pdf":
        with PdfPages(str(output_path)) as pdf:
            pdf.savefig(fig, bbox_inches="tight", dpi=300)
    else:
        fig.savefig(str(output_path), format=fmt, bbox_inches="tight", dpi=300)


def render_logo_plot(
    title: str,
    matrix: pd.DataFrame,
    labels: list[str],
    output_path: Path,
    color_scheme: str = "chemistry",
    matrix_type: str = "information",
    region_groups: dict[str, list[int]] | None = None,
) -> None:
    """Render logo plot(s) and save.

    Single panel for <=15 positions or no region grouping.
    Multi-panel (one subplot per region) for >15 positions with regions.
    """
    apply_nature_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    y_label = _y_axis_label(matrix_type)
    n_positions = len(labels)
    fmt = output_path.suffix.lstrip(".") or "pdf"

    use_multi = (
        region_groups is not None
        and len(region_groups) > 1
        and n_positions > 15
    )

    if not use_multi:
        fig_width = min(
            max(NATURE_SINGLE_COL, n_positions * 0.5), NATURE_DOUBLE_COL
        )
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, PANEL_HEIGHT))

        logo = logomaker.Logo(
            matrix,
            ax=ax,
            color_scheme=color_scheme,
            font_name="Arial",
            stack_order="big_on_top",
            vpad=0.04,
            baseline_width=0.4,
        )
        logo.style_spines(spines=["top", "right"], visible=False)
        logo.style_spines(spines=["left", "bottom"], visible=True)

        ax.set_xticks(range(n_positions))
        ax.set_xticklabels(labels, rotation=90, fontsize=FONT_TICK_LABEL)
        ax.set_ylabel(y_label, fontsize=FONT_AXIS_LABEL, labelpad=4)
        ax.set_title(
            f"{title} \u2014 {n_positions} selected positions",
            fontsize=FONT_TITLE, fontweight="bold", pad=4, loc="left",
        )
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=3, prune="both")
        )

        plt.tight_layout(pad=0.4)
        _save_figure(fig, output_path, fmt)
        plt.close(fig)

    else:
        sorted_regions = [
            name for name, _, _ in ENV_REGIONS if name in region_groups
        ]
        n_regions = len(sorted_regions)
        fig, axes = plt.subplots(
            n_regions, 1,
            figsize=(
                NATURE_DOUBLE_COL,
                max(3.0, PANEL_HEIGHT_MULTI * n_regions),
            ),
            squeeze=False,
        )

        for idx, region_name in enumerate(sorted_regions):
            ax = axes[idx, 0]
            row_indices = region_groups[region_name]
            sub_labels = [labels[r] for r in row_indices]
            sub_matrix = matrix.iloc[row_indices, :].copy()
            sub_matrix.index = range(len(row_indices))

            logo = logomaker.Logo(
                sub_matrix,
                ax=ax,
                color_scheme=color_scheme,
                font_name="Arial",
                stack_order="big_on_top",
                vpad=0.04,
                baseline_width=0.4,
            )
            logo.style_spines(spines=["top", "right"], visible=False)
            logo.style_spines(spines=["left", "bottom"], visible=True)

            ax.set_xticks(range(len(row_indices)))
            ax.set_xticklabels(
                sub_labels, rotation=90, fontsize=FONT_TICK_LABEL
            )
            ax.set_ylabel(y_label, fontsize=FONT_AXIS_LABEL, labelpad=4)
            ax.set_title(
                f"{title} \u2014 {region_name}",
                fontsize=FONT_TITLE, fontweight="bold", pad=4, loc="left",
            )
            ax.yaxis.set_major_locator(
                ticker.MaxNLocator(nbins=3, prune="both")
            )

        plt.tight_layout(pad=0.4, h_pad=0.6)
        _save_figure(fig, output_path, fmt)
        plt.close(fig)


def process_standard_logo(
    alignment: dict[str, str],
    hxb2_map: list[HxB2Position],
    sample_seqs: list[str],
    output_path: Path,
    title: str = "",
    positions: list[int] | None = None,
    regions: list[str] | None = None,
    matrix_type: str = "information",
    color_scheme: str = "chemistry",
    max_gap_fraction: float = 0.9,
) -> bool:
    """Orchestrator for standard logo mode.

    Resolves positions, filters gaps, builds matrix, renders.

    Returns:
        True if output was produced.
    """
    reverse_map = build_reverse_hxb2_map(hxb2_map)

    # Resolve HxB2 positions using priority order
    hxb2_positions: list[int] = []
    active_regions: list[str] | None = None

    if positions is not None:
        hxb2_positions = [p for p in positions if p in reverse_map]
    elif regions is not None:
        hxb2_positions = resolve_positions_from_regions(regions, reverse_map)
        active_regions = regions
    else:
        hxb2_positions = resolve_positions_from_regions(
            V_LOOP_REGIONS, reverse_map
        )
        active_regions = V_LOOP_REGIONS

    if not hxb2_positions:
        return False

    # Filter by gap fraction
    selected_cols, valid_hxb2 = filter_positions_by_gaps(
        hxb2_positions, reverse_map, sample_seqs, max_gap_fraction
    )

    if not selected_cols:
        return False

    # Build region groups for multi-panel
    region_groups = build_region_groups(valid_hxb2)

    # Build logo matrix
    matrix, labels = build_logo_matrix(
        sample_seqs, selected_cols, hxb2_map, matrix_type
    )

    # Render
    render_logo_plot(
        title,
        matrix,
        labels,
        output_path,
        color_scheme,
        matrix_type,
        region_groups=region_groups if len(region_groups) > 1 else None,
    )

    return True
