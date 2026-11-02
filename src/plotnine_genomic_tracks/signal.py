from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from plotnine import geom_path, geom_rect
from pygenomic_tracks import as_region
from pygenomic_tracks._coords import DEFAULT_COLS
from pygenomic_tracks.geometry import signal_steps

from ._common import mapped, mapped_columns, with_mapping

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from plotnine import aes
    from pygenomic_tracks._coords import RegionLike

__all__ = ["geom_signal", "signal_data"]


def signal_data(
    data: pd.DataFrame,
    region: RegionLike,
    value: str = "signal",
    *,
    by: Sequence[str] | None = None,
    cols: Sequence[str] = DEFAULT_COLS,
    n_bins: int | None = 2000,
    agg: str = "max",
) -> pd.DataFrame:
    """
    Steps of a signal track, one table per group stacked together

    The computation behind `geom_signal`, for inspecting or plotting
    yourself.

    Returns
    -------
    pandas.DataFrame
        ``start``, ``end``, ``value``, plus the ``by`` columns.
    """
    region = as_region(region)
    by = list(by or [])
    parts = []
    groups = data.groupby(by, sort=False, observed=True) if by else [((), data)]
    for key, sub in groups:
        steps = signal_steps(
            sub, region, value=value, cols=cols, n_bins=n_bins, agg=agg
        )
        key = key if isinstance(key, tuple) else (key,)
        for col, val in zip(by, key):
            steps[col] = val
        parts.append(steps)
    if not parts:
        return pd.DataFrame(columns=["start", "end", "value", *by])
    return pd.concat(parts, ignore_index=True)


def _stair_paths(steps: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """
    Vertices of the outline of each contiguous run of steps
    """
    rows = []
    groups = steps.groupby(by, sort=False, observed=True) if by else [((), steps)]
    for k, (_, g) in enumerate(groups):
        s = g["start"].to_numpy(float)
        e = g["end"].to_numpy(float)
        v = g["value"].to_numpy(float)
        run = np.concatenate([[0], np.cumsum(s[1:] != e[:-1])])
        frame = pd.DataFrame(
            {
                "x": np.column_stack([s, e]).ravel(),
                "y": np.repeat(v, 2),
                "_path": [f"{k}:{r}" for r in np.repeat(run, 2)],
            }
        )
        for col in by:
            frame[col] = g[col].iloc[0]
        rows.append(frame)
    return pd.concat(rows, ignore_index=True)


def geom_signal(
    data: pd.DataFrame,
    region: RegionLike,
    value: str = "signal",
    mapping: aes | None = None,
    *,
    by: Sequence[str] | None = None,
    cols: Sequence[str] = DEFAULT_COLS,
    n_bins: int | None = 2000,
    agg: str = "max",
    fill: bool = True,
    color: str = "#3b6ea8",
    size: float = 0.4,
    **kwargs: Any,
):
    """
    Quantitative signal track: coverage, ChEC-seq, ChIP-seq, bigWig values

    Each interval is drawn as a step over its real extent; gaps and NaN
    values stay empty. With more intervals than ``n_bins``, or overlapping
    intervals, they are aggregated onto ``n_bins`` bins across the region
    (``"max"`` by default, so narrow peaks survive). The computation is
    `pygenomic_tracks.geometry.signal_steps`, the same as in
    pygenomic_tracks' matplotlib tracks.

    Parameters
    ----------
    data : pandas.DataFrame
        Intervals with ``chrom``, ``start``, ``end`` (0-based, half-open)
        and a numeric ``value`` column, such as ``pymisha.gextract``
        output.
    region : Region, str or tuple
        The window, e.g. ``"chr2:1,200,000-1,350,000"``.
    value : str
        Column holding the signal.
    mapping : aes, optional
        Extra aesthetics, e.g. ``aes(fill="sample")`` to overlay samples.
        Columns it names are binned separately.
    by : sequence of str, optional
        More columns to bin separately and keep, e.g. a column you facet
        on.
    cols : tuple of str
        Names of the chrom, start and end columns.
    n_bins : int or None
        Resolution of the aggregation. 2000 is finer than any panel is
        wide in pixels. `None` never aggregates and raises on
        overlapping intervals.
    agg : {"max", "min", "mean"}
        Aggregation per bin; ``"mean"`` is weighted by overlap.
    fill : bool
        Filled area (`True`) or outline.
    color : str
        Fill or line colour, unless mapped.
    size : float
        Outline width when ``fill=False``.
    **kwargs
        Passed to `geom_rect` (filled) or `geom_path` (outline), e.g.
        ``alpha``.

    Returns
    -------
    geom
        A plotnine layer; add it with ``+``.
    """
    by = list(by or [])
    by += [c for c in mapped_columns(mapping, data.columns) if c not in by]
    steps = signal_data(
        data, region, value, by=by, cols=cols, n_bins=n_bins, agg=agg
    )
    user = mapped(mapping)

    if fill:
        params = {} if {"fill"} & user else {"fill": color}
        return geom_rect(
            with_mapping(
                {"xmin": "start", "xmax": "end", "ymin": 0, "ymax": "value"},
                mapping,
            ),
            data=steps,
            inherit_aes=False,
            **params,
            **kwargs,
        )

    params = {} if {"color", "colour"} & user else {"color": color}
    return geom_path(
        with_mapping({"x": "x", "y": "y", "group": "_path"}, mapping),
        data=_stair_paths(steps, by),
        inherit_aes=False,
        size=size,
        **params,
        **kwargs,
    )
