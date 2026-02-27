"""Shared constants for tlogo: Env regions, Nature typography, color schemes."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Env gp160 region boundaries (HxB2 AA numbering, LANL convention)
# ---------------------------------------------------------------------------

ENV_REGIONS: list[tuple[str, int, int]] = [
    ("SP",   1,  30),    # Signal peptide
    ("C1",  31, 130),
    ("V1", 131, 157),    # C131-C157 disulfide bounded
    ("V2", 158, 196),    # S158-C196
    ("C2", 197, 295),
    ("V3", 296, 331),    # C296-C331 disulfide bounded
    ("C3", 332, 384),
    ("V4", 385, 418),    # C385-C418 disulfide bounded
    ("C4", 419, 459),
    ("V5", 460, 469),    # N460-R469
    ("C5", 470, 511),    # includes REKR cleavage site at 508-511
    ("gp41", 512, 856),
]

V_LOOP_REGIONS: list[str] = ["V1", "V2", "V3", "V4", "V5"]

MAX_GAP_FRACTION: float = 0.9

COLOR_SCHEMES: list[str] = ["chemistry", "charge", "hydrophobicity"]

# ---------------------------------------------------------------------------
# Nature figure typography constants (pt)
# ---------------------------------------------------------------------------

FONT_TICK_LABEL: int = 6    # X-axis HxB2 position labels
FONT_AXIS_LABEL: int = 7    # Y-axis "Information (bits)" label
FONT_TITLE: int = 8         # Panel title (single panel or per-region title)
FONT_SUPTITLE: int = 8      # Figure-level suptitle in window mode
FONT_ANIMAL_LABEL: int = 7  # Per-animal label in window mode

# ---------------------------------------------------------------------------
# Nature figure sizing constants (inches)
# ---------------------------------------------------------------------------

NATURE_SINGLE_COL: float = 3.5   # 89 mm -- Nature single column
NATURE_DOUBLE_COL: float = 7.2   # 183 mm -- Nature double column
PANEL_HEIGHT: float = 1.8         # Height per logo panel (single or window)
PANEL_HEIGHT_MULTI: float = 1.6   # Height per region panel in multi-panel layout
