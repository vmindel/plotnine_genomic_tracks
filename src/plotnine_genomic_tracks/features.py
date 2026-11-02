from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from plotnine import geom_point, geom_rect, geom_segment, geom_text
from pygenomic_tracks import annotation, as_region
from pygenomic_tracks._coords import DEFAULT_COLS, check_columns
from pygenomic_tracks.geometry import (
    clipped,
    feature_rows,
    gene_models,
    select,
    widen,
)

from ._common import (
    LABEL_PITCH,
    PLAIN_PITCH,
    bp_per_px,
    chevrons,
    feature_axis,
    mapped,
    row_limits,
    with_mapping,
)

if TYPE_CHECKING:
    import os
    from collections.abc import Sequence
    from typing import Any

    from plotnine import aes
    from pygenomic_tracks._coords import RegionLike

__all__ = ["gene_data", "geom_genes", "geom_intervals", "interval_data"]

# Strand arrowheads: one per this many pixels, none on features narrower
# than MIN_CHEVRON_PX. Matches pygenomic_tracks.
CHEVRON_SPACING_PX = 40.0
MIN_CHEVRON_PX = 6.0


def _chevron_layers(xs, ys, strands, color, size):
    layers = []
    for st, shape in (("+", ">"), ("-", "<")):
        keep = strands == st
        if keep.any():
            layers.append(
                geom_point(
                    with_mapping({"x": "x", "y": "y"}, None),
                    data=pd.DataFrame({"x": xs[keep], "y": ys[keep]}),
                    inherit_aes=False,
                    shape=shape,
                    color=color,
                    fill=color,
                    size=size,
                    stroke=0,
                )
            )
    return layers


def interval_data(
    data: pd.DataFrame,
    region: RegionLike,
    *,
    cols: Sequence[str] = DEFAULT_COLS,
    label_by: str | None = None,
    pack: bool = True,
    height: float = 0.6,
    min_width_px: float = 1.0,
    fontsize: float = 7,
    panel_width: float = 6.0,
    dpi: float = 100,
) -> pd.DataFrame:
    """
    Rows of `data` in the region, with the box and row they are drawn in

    The computation behind `geom_intervals`. Adds ``xmin``, ``xmax``
    (clipped, widened to at least ``min_width_px``), ``row``, ``y``,
    ``ymin``, ``ymax`` and ``xmid``.
    """
    region = as_region(region)
    if label_by is not None:
        check_columns(data, [label_by])
    sub = select(data, region, cols).copy()
    bpp = bp_per_px(region, panel_width, dpi)
    starts, ends = clipped(sub, region, cols)
    starts, ends = widen(starts, ends, min_width_px * bpp)
    labels = sub[label_by].astype(str).tolist() if label_by is not None else None
    rows = feature_rows(
        starts, ends, bpp, labels=labels, fontsize=fontsize, dpi=dpi, pack=pack
    )
    pitch = LABEL_PITCH if label_by is not None else PLAIN_PITCH
    y = -rows * pitch
    sub["xmin"], sub["xmax"] = starts, ends
    sub["xmid"] = (starts + ends) / 2
    sub["row"] = rows
    sub["y"] = y
    sub["ymin"], sub["ymax"] = y - height / 2, y + height / 2
    return sub.reset_index(drop=True)


def geom_intervals(
    data: pd.DataFrame,
    region: RegionLike,
    mapping: aes | None = None,
    *,
    cols: Sequence[str] = DEFAULT_COLS,
    label_by: str | None = None,
    strand: str | None = None,
    pack: bool = True,
    height: float = 0.6,
    min_width_px: float = 1.0,
    fill: str = "#595959",
    color: str | None = None,
    fontsize: float = 7,
    panel_width: float = 6.0,
    dpi: float = 100,
    **kwargs: Any,
) -> list:
    """
    Feature track of boxes: peaks, CREs, motif hits, any BED-like table

    Overlapping features, and their labels, are packed into rows. Colour
    is an ordinary plotnine aesthetic: ``aes(fill="cre")`` with any
    ``scale_fill_*`` gives a legend.

    Parameters
    ----------
    data : pandas.DataFrame
        Intervals with ``chrom``, ``start``, ``end``, 0-based half-open.
    region : Region, str or tuple
        The window to draw.
    mapping : aes, optional
        Aesthetics for the boxes, e.g. ``aes(fill="cre")`` or
        ``aes(fill="score")``.
    cols : tuple of str
        Names of the chrom, start and end columns.
    label_by : str, optional
        Column printed under each feature; packing leaves room for it.
    strand : str, optional
        Strand column (``"+"``/``"-"``); draws arrowheads on features at
        least 6 px wide.
    pack : bool
        Stack overlapping features into rows; `False` draws one row.
    height : float
        Box height as a fraction of a row.
    min_width_px : float
        Features narrower than this are widened around their centre, so
        10 bp motifs stay visible in a 1 Mb window.
    fill, color : str
        Box fill and outline, unless mapped.
    fontsize : float
        Label size in points.
    panel_width, dpi : float
        Panel width (inches) and resolution assumed for packing, minimum
        widths and arrowhead spacing; see the docs on resolution.
    **kwargs
        Passed to the box layer (`geom_rect`), e.g. ``alpha``.

    Returns
    -------
    list
        Layers plus the track's y scale and theme; add with ``+``.
    """
    region = as_region(region)
    if strand is not None:
        check_columns(data, [strand])
    frame = interval_data(
        data,
        region,
        cols=cols,
        label_by=label_by,
        pack=pack,
        height=height,
        min_width_px=min_width_px,
        fontsize=fontsize,
        panel_width=panel_width,
        dpi=dpi,
    )
    params = {} if {"fill"} & mapped(mapping) else {"fill": fill}
    # Pass a colour only when set: an explicit None breaks faceted plots
    # in some plotnine versions.
    if color is not None and not {"color", "colour"} & mapped(mapping):
        params["color"] = color
    layers: list = [
        geom_rect(
            with_mapping(
                {"xmin": "xmin", "xmax": "xmax", "ymin": "ymin", "ymax": "ymax"},
                mapping,
            ),
            data=frame,
            inherit_aes=False,
            **params,
            **kwargs,
        )
    ]
    if strand is not None and len(frame):
        bpp = bp_per_px(region, panel_width, dpi)
        xs, ys, ss = chevrons(
            frame["xmin"],
            frame["xmax"],
            frame["y"],
            frame[strand].astype(str),
            CHEVRON_SPACING_PX * bpp,
            MIN_CHEVRON_PX * bpp,
        )
        layers += _chevron_layers(xs, ys, ss, "white", size=1.6)
    if label_by is not None and len(frame):
        layers.append(
            geom_text(
                with_mapping({"x": "xmid", "y": "_label_y", "label": label_by}, None),
                data=frame.assign(_label_y=frame["ymin"] - 0.08),
                inherit_aes=False,
                va="top",
                size=fontsize,
                color="#262626",
            )
        )
    n_rows = int(frame["row"].max()) + 1 if len(frame) else 1
    pitch = LABEL_PITCH if label_by is not None else PLAIN_PITCH
    return layers + feature_axis(
        row_limits(n_rows, pitch, height, label_by is not None)
    )


def gene_data(
    genes: pd.DataFrame | str | os.PathLike,
    region: RegionLike,
    *,
    cols: Sequence[str] = DEFAULT_COLS,
    feature: str = "feature",
    transcript: str = "transcript_id",
    name: str | None = "gene_name",
    strand: str | None = "strand",
    transcripts: str = "all",
    height: float = 0.6,
    utr_ratio: float = 0.5,
    show_names: bool = True,
    pack: bool = True,
    fontsize: float = 7,
    panel_width: float = 6.0,
    dpi: float = 100,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Transcripts and their parts, positioned in rows

    The computation behind `geom_genes`.

    Returns
    -------
    spans : pandas.DataFrame
        One row per transcript: ``start``, ``end`` (clipped), ``name``,
        ``strand``, ``row``, ``y``.
    parts : pandas.DataFrame
        Exon/CDS boxes: ``transcript``, ``xmin``, ``xmax``, ``ymin``,
        ``ymax``, ``feature``, ``full``.
    """
    region = as_region(region)
    df = annotation(genes)
    span, parts = gene_models(
        df,
        region,
        cols=cols,
        feature=feature,
        transcript=transcript,
        name=name,
        strand=strand,
        transcripts=transcripts,
    )
    span = span.copy()
    span["start"] = np.clip(span["start"].to_numpy(float), region.start, region.end)
    span["end"] = np.clip(span["end"].to_numpy(float), region.start, region.end)
    bpp = bp_per_px(region, panel_width, dpi)
    labels = span["name"].astype(str).tolist() if show_names else None
    rows = feature_rows(
        span["start"], span["end"], bpp, labels=labels, fontsize=fontsize,
        dpi=dpi, pack=pack,
    )
    pitch = LABEL_PITCH if show_names else PLAIN_PITCH
    span["row"] = rows
    span["y"] = -rows * pitch

    parts = parts.copy()
    y_of = span["y"].to_dict()
    h = np.where(parts["full"], height, height * utr_ratio)
    yc = parts["transcript"].map(y_of).to_numpy(float)
    parts["xmin"] = np.clip(parts["start"].to_numpy(float), region.start, region.end)
    parts["xmax"] = np.clip(parts["end"].to_numpy(float), region.start, region.end)
    parts["ymin"], parts["ymax"] = yc - h / 2, yc + h / 2
    parts = parts[parts["xmax"] > parts["xmin"]]
    span.index.name = "transcript"
    return span.reset_index(), parts.reset_index(drop=True)


def geom_genes(
    genes: pd.DataFrame | str | os.PathLike,
    region: RegionLike,
    *,
    transcripts: str = "all",
    cols: Sequence[str] = DEFAULT_COLS,
    feature: str = "feature",
    transcript: str = "transcript_id",
    name: str | None = "gene_name",
    strand: str | None = "strand",
    color: str = "#333333",
    height: float = 0.6,
    utr_ratio: float = 0.5,
    show_names: bool = True,
    pack: bool = True,
    fontsize: float = 7,
    panel_width: float = 6.0,
    dpi: float = 100,
) -> list:
    """
    Gene model track: intron line, exons, thicker CDS, strand arrowheads, names

    Parameters
    ----------
    genes : pandas.DataFrame, str or path
        A UCSC assembly name (``"hg38"``, ``"mm10"``; downloaded and cached
        by pygenomic_tracks on first use), a GTF path, or a DataFrame with
        one row per exon/CDS such as `pygenomic_tracks.read_gtf` returns.
    region : Region, str or tuple
        The window to draw. Transcripts crossing its edges are clipped.
    transcripts : {"all", "longest"}
        ``"longest"`` keeps one transcript per gene name.
    cols, feature, transcript, name, strand : str
        Column names, as in `pygenomic_tracks.genes`. ``name=None`` shows
        transcript ids; ``strand=None`` draws no arrowheads.
    color : str
        Colour of the whole model.
    height : float
        Height of coding exons as a fraction of a row.
    utr_ratio : float
        Height of UTRs relative to ``height``.
    show_names : bool
        Print gene names under transcripts.
    pack : bool
        Stack overlapping transcripts into rows.
    fontsize : float
        Size of the names, in points.
    panel_width, dpi : float
        Panel width (inches) and resolution assumed for packing and
        arrowhead spacing.

    Returns
    -------
    list
        Layers plus the track's y scale and theme; add with ``+``.
    """
    region = as_region(region)
    spans, parts = gene_data(
        genes,
        region,
        cols=cols,
        feature=feature,
        transcript=transcript,
        name=name,
        strand=strand,
        transcripts=transcripts,
        height=height,
        utr_ratio=utr_ratio,
        show_names=show_names,
        pack=pack,
        fontsize=fontsize,
        panel_width=panel_width,
        dpi=dpi,
    )
    pitch = LABEL_PITCH if show_names else PLAIN_PITCH
    n_rows = int(spans["row"].max()) + 1 if len(spans) else 1
    axis = feature_axis(row_limits(n_rows, pitch, height, show_names))
    if not len(spans):
        return axis

    layers: list = [
        geom_segment(
            with_mapping({"x": "start", "xend": "end", "y": "y", "yend": "y"}, None),
            data=spans,
            inherit_aes=False,
            color=color,
            size=0.45,
        ),
        geom_rect(
            with_mapping(
                {"xmin": "xmin", "xmax": "xmax", "ymin": "ymin", "ymax": "ymax"},
                None,
            ),
            data=parts,
            inherit_aes=False,
            fill=color,
        ),
    ]
    if strand is not None:
        bpp = bp_per_px(region, panel_width, dpi)
        xs, ys, ss = chevrons(
            spans["start"],
            spans["end"],
            spans["y"],
            spans["strand"].astype(str),
            CHEVRON_SPACING_PX * bpp,
            MIN_CHEVRON_PX * bpp,
        )
        layers += _chevron_layers(xs, ys, ss, color, size=1.8)
    if show_names:
        layers.append(
            geom_text(
                with_mapping({"x": "_mid", "y": "_label_y", "label": "name"}, None),
                data=spans.assign(
                    _mid=(spans["start"] + spans["end"]) / 2,
                    _label_y=spans["y"] - height / 2 - 0.08,
                ),
                inherit_aes=False,
                va="top",
                size=fontsize,
                color="#262626",
            )
        )
    return layers + axis
