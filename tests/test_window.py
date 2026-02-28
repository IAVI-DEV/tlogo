"""Tests for tlogo.window module (rendering + orchestration)."""

from __future__ import annotations

from pathlib import Path

from tlogo.hxb2 import HxB2Position, build_hxb2_map, build_reverse_hxb2_map
from tlogo.window import process_window_mode, render_window_logo


class TestRenderWindowLogo:
    def test_renders_without_error(self, tmp_path):
        """render_window_logo produces a file."""
        hxb2_map = [HxB2Position(i, i + 1, "SP", "A") for i in range(10)]
        animal_groups = {
            "animalA": ["AAAAAAAAAA", "AAAAABAAAA"],
            "animalB": ["AAAAAAAAAA", "AAAAACAAAA"],
        }
        selected_cols = [3, 4, 5, 6, 7]
        labels = ["4", "5", "6", "7", "8"]
        out = tmp_path / "test_window.pdf"

        render_window_logo(
            "Test", 5, animal_groups, selected_cols, labels,
            hxb2_map, out, matrix_type="information",
        )
        assert out.exists()
        assert out.stat().st_size > 0


class TestProcessWindowMode:
    def _make_alignment(self):
        seq = "MRVKEKYQHL"
        return {
            "HxB2": seq,
            "parental_ref": seq,
            "animalA_wk8_001": seq,
            "animalA_wk8_002": "MRVKEKYQRL",
        }

    def test_returns_true_on_success(self, tmp_path):
        alignment = self._make_alignment()
        hxb2_map = build_hxb2_map(alignment)
        animal_groups = {
            "animalA": [alignment["animalA_wk8_001"], alignment["animalA_wk8_002"]],
        }
        out = tmp_path / "window.pdf"

        result = process_window_mode(
            alignment, hxb2_map, 5, 2, animal_groups, out,
            title="test",
        )
        assert result is True
        assert out.exists()

    def test_returns_false_for_missing_position(self, tmp_path):
        alignment = self._make_alignment()
        hxb2_map = build_hxb2_map(alignment)
        animal_groups = {
            "animalA": [alignment["animalA_wk8_001"]],
        }
        out = tmp_path / "window.pdf"

        result = process_window_mode(
            alignment, hxb2_map, 9999, 2, animal_groups, out,
        )
        assert result is False

    def test_returns_false_for_empty_groups(self, tmp_path):
        alignment = self._make_alignment()
        hxb2_map = build_hxb2_map(alignment)
        out = tmp_path / "window.pdf"

        result = process_window_mode(
            alignment, hxb2_map, 5, 2, {}, out,
        )
        assert result is False
