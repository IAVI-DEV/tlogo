"""Window-mode and auto-variants logo rendering."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import logomaker

from tlogo.constants import (
    FONT_ANIMAL_LABEL,
    FONT_SUPTITLE,
    NATURE_DOUBLE_COL,
    NATURE_SINGLE_COL,
    PANEL_HEIGHT,
)
from tlogo.hxb2 import HxB2Position, build_reverse_hxb2_map, get_env_region
from tlogo.io import sort_animal_groups
from tlogo.logo import (
    apply_nature_style,
    draw_logo_on_axes,
    save_figure,
    y_axis_label,
)
from tlogo.matrix import find_variant_positions, resolve_window_cols


def render_window_logo(
    title: str,
    center_hxb2: int,
    animal_groups: dict[str, list[str]],
    selected_cols: list[int],
    labels: list[str],
    hxb2_map: list[HxB2Position],
    output_path: Path,
    color_scheme: str = "chemistry",
    matrix_type: str = "information",
    self_name: str | None = None,
) -> None:
    """Render per-animal logo subplots for a window around one HxB2 position.

    One subplot row per animal, stacked vertically. Each logo shows the
    amino acid frequency distribution for that animal's sequences only.
    """
    apply_nature_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    y_label = y_axis_label(matrix_type)
    pseudo = 1.0 if matrix_type == "information" else 0.0
    fmt = output_path.suffix.lstrip(".") or "pdf"

    sorted_animals = sort_animal_groups(
        list(animal_groups.keys()), self_name
    )
    n_animals = len(sorted_animals)

    fig_width = min(
        max(NATURE_SINGLE_COL, len(selected_cols) * 0.8), NATURE_DOUBLE_COL
    )
    fig_height = max(3, PANEL_HEIGHT * n_animals)
    fig, axes = plt.subplots(
        n_animals, 1, figsize=(fig_width, fig_height), squeeze=False,
    )

    region = get_env_region(center_hxb2) or ""
    fig.suptitle(
        f"{title} \u2014 HxB2 {center_hxb2} ({region}) \u00b1 window",
        fontsize=FONT_SUPTITLE, fontweight="bold", y=0.98,
    )

    for idx, animal_name in enumerate(sorted_animals):
        ax = axes[idx, 0]
        seqs = animal_groups[animal_name]
        n_seqs = len(seqs)

        subsequences = [
            "".join(s[c] for c in selected_cols) for s in seqs
        ]

        matrix = logomaker.alignment_to_matrix(
            subsequences,
            to_type=matrix_type,
            characters_to_ignore=".-",
            pseudocount=pseudo,
        )
        matrix.index = range(len(matrix))

        draw_logo_on_axes(
            ax, matrix, labels,
            title_text=f"{animal_name}  (n={n_seqs})",
            y_label=y_label,
            color_scheme=color_scheme,
            title_fontsize=FONT_ANIMAL_LABEL,
        )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, output_path, fmt)
    plt.close(fig)


def process_window_mode(
    alignment: dict[str, str],
    hxb2_map: list[HxB2Position],
    window_pos: int,
    radius: int,
    animal_groups: dict[str, list[str]],
    output_path: Path,
    title: str = "",
    matrix_type: str = "information",
    color_scheme: str = "chemistry",
    self_name: str | None = None,
) -> bool:
    """Orchestrator for window mode on a single alignment.

    Returns:
        True if a PDF was generated.
    """
    reverse_map = build_reverse_hxb2_map(hxb2_map)
    all_sample_seqs = [s for seqs in animal_groups.values() for s in seqs]

    if not all_sample_seqs:
        return False

    if window_pos not in reverse_map:
        return False

    selected_cols, labels = resolve_window_cols(
        window_pos, radius, reverse_map, hxb2_map, all_sample_seqs,
    )

    if not selected_cols:
        return False

    render_window_logo(
        title, window_pos, animal_groups, selected_cols, labels,
        hxb2_map, output_path, color_scheme, matrix_type, self_name,
    )
    return True


def process_auto_variants(
    alignment: dict[str, str],
    hxb2_map: list[HxB2Position],
    animal_groups: dict[str, list[str]],
    ref_id: str,
    output_dir: Path,
    title: str = "",
    radius: int = 4,
    variant_freq: float = 0.50,
    matrix_type: str = "information",
    color_scheme: str = "chemistry",
    self_name: str | None = None,
    fmt: str = "pdf",
) -> list[int]:
    """Orchestrator for auto-variants mode.

    Returns:
        List of HxB2 positions where window logos were generated.
    """
    var_positions = find_variant_positions(
        alignment, ref_id, hxb2_map, animal_groups, variant_freq,
    )

    generated: list[int] = []
    for vpos in var_positions:
        out_path = output_dir / f"logo_window_{vpos}.{fmt}"
        if process_window_mode(
            alignment, hxb2_map, vpos, radius, animal_groups,
            out_path, title, matrix_type, color_scheme, self_name,
        ):
            generated.append(vpos)

    return generated
