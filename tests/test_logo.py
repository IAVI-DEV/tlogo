"""Tests for tlogo.logo module (rendering helpers)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from tlogo.logo import (
    apply_nature_style,
    draw_logo_on_axes,
    save_figure,
    y_axis_label,
)


class TestYAxisLabel:
    def test_information(self):
        assert y_axis_label("information") == "Information (bits)"

    def test_probability(self):
        assert y_axis_label("probability") == "Frequency"

    def test_counts(self):
        assert y_axis_label("counts") == "Count"

    def test_unknown(self):
        assert y_axis_label("other") == "Value"


class TestApplyNatureStyle:
    def test_sets_agg_backend(self):
        apply_nature_style()
        import matplotlib
        assert matplotlib.get_backend().lower() == "agg"


class TestSaveFigure:
    def test_save_pdf(self, tmp_path):
        apply_nature_style()
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        out = tmp_path / "test.pdf"
        save_figure(fig, out, "pdf")
        plt.close(fig)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_save_png(self, tmp_path):
        apply_nature_style()
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        out = tmp_path / "test.png"
        save_figure(fig, out, "png")
        plt.close(fig)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_creates_parent_dirs(self, tmp_path):
        apply_nature_style()
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        out = tmp_path / "nested" / "dir" / "test.pdf"
        save_figure(fig, out, "pdf")
        plt.close(fig)
        assert out.exists()


class TestDrawLogoOnAxes:
    def test_renders_without_error(self, tmp_path):
        """draw_logo_on_axes produces a styled axes."""
        import logomaker

        apply_nature_style()
        seqs = ["ACDE", "ACDE", "ADDE"]
        matrix = logomaker.alignment_to_matrix(
            seqs, to_type="information", characters_to_ignore=".-",
        )
        matrix.index = range(len(matrix))

        fig, ax = plt.subplots()
        draw_logo_on_axes(
            ax, matrix, ["1", "2", "3", "4"],
            title_text="Test title",
            y_label="Information (bits)",
        )

        # Title is set with loc="left", so check _left_title
        assert ax.get_title(loc="left") == "Test title"
        assert ax.get_ylabel() == "Information (bits)"
        plt.close(fig)
