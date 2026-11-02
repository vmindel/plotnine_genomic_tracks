from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from plotnine import aes, element_blank, scale_y_continuous, theme

if TYPE_CHECKING:
    from typing import Any

    from pygenomic_tracks import Region

# Row pitch, in units of box height, as in pygenomic_tracks' matplotlib tracks
LABEL_PITCH = 1.8
PLAIN_PITCH = 1.0


def bp_per_px(region: Region, panel_width: float, dpi: float) -> float:
    """
    Resolution assumed when packing and spacing features

    plotnine only knows the panel size when drawing, after the rows of a
    feature track must already be fixed. Packing therefore assumes a panel
    ``panel_width`` inches wide at ``dpi``.
    """
    return region.length / max(1.0, panel_width * dpi)


def with_mapping(ours: dict[str, Any], user: aes | None) -> aes:
    """
    Our aesthetics plus the user's, the user's taking precedence
    """
    return aes(**{**ours, **dict(user or {})})


def mapped(user: aes | None) -> set[str]:
    return set(dict(user or {}))


def mapped_columns(user: aes | None, columns) -> list[str]:
    """
    Data columns named by the user's mapping, in mapping order
    """
    out = []
    for value in dict(user or {}).values():
        if isinstance(value, str) and value in columns and value not in out:
            out.append(value)
    return out


def chevrons(starts, ends, ys, strands, step_bp: float, min_bp: float):
    """
    Positions of strand arrowheads along features

    Returns x, y and strand ("+" or "-") of each arrowhead. Features shorter
    than ``min_bp`` and features without a strand get none.
    """
    xs, yy, ss = [], [], []
    for s, e, y, st in zip(starts, ends, ys, strands):
        if st not in ("+", "-") or (e - s) < min_bp:
            continue
        n = max(1, int((e - s) // step_bp))
        xs.extend(s + (np.arange(n) + 0.5) * (e - s) / n)
        yy.extend([y] * n)
        ss.extend([st] * n)
    return np.asarray(xs, float), np.asarray(yy, float), np.asarray(ss, object)


def row_limits(n_rows: int, pitch: float, height: float, labelled: bool):
    """
    y-limits for feature rows drawn at y = -row * pitch
    """
    bottom = (max(n_rows, 1) - 1) * pitch + height / 2 + (
        pitch - height if labelled else 0.2
    )
    return (-bottom, height / 2 + 0.2)


def feature_axis(limits) -> list:
    """
    Scale and theme of a feature track: fixed rows, no y-axis
    """
    return [
        scale_y_continuous(limits=limits, breaks=[], expand=(0, 0)),
        theme(
            axis_line_y=element_blank(),
            axis_ticks_y=element_blank(),
            axis_text_y=element_blank(),
        ),
    ]
