import io

import numpy as np
import pandas as pd
import pygenomic_tracks as pygt
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from plotnine import aes, facet_wrap, geom_path, geom_rect, geom_text, ggplot, labs
from plotnine.exceptions import PlotnineError
from pygenomic_tracks import geometry as G

import plotnine_genomic_tracks as pgt

REGION = pygt.example_region()


def render(p):
    buf = io.BytesIO()
    p.save(buf, format="png", width=6, height=3, dpi=72, verbose=False)
    return buf.getbuffer().nbytes




# scale_x_genome


def test_scale_x_genome_units_and_name():
    scale, coord = pgt.scale_x_genome(REGION)
    assert scale.name == "chr2 (kb)"
    assert scale.labels([1_200_000, 1_250_500]) == ["1,200", "1,250.5"]
    assert coord.limits.x == (REGION.start, REGION.end)
    scale, _ = pgt.scale_x_genome(REGION, unit="Mb", name=False)
    assert scale.labels([1_250_000]) == ["1.25"] and scale.name is None
    with pytest.raises(ValueError):
        pgt.scale_x_genome(REGION, unit="cM")


def test_scale_x_genome_keeps_data_outside_region():
    d = pd.DataFrame({"x": [REGION.start - 10_000, REGION.end + 10_000], "y": [0, 1]})
    p = ggplot(d) + geom_path(aes("x", "y")) + pgt.scale_x_genome(REGION)
    fig = p.draw()
    FigureCanvasAgg(fig).draw()
    ax = fig.axes[0]
    assert ax.get_xlim() == (REGION.start, REGION.end)
    drawn = [c for c in ax.collections if len(c.get_paths())] + list(ax.lines)
    points = sum(len(getattr(a, "get_xydata", lambda: a.get_paths()[0].vertices)()) for a in drawn)
    assert points == 2


# Signal


def test_signal_matches_pygenomic_tracks_geometry():
    sig = pygt.example_signal()
    layer = pgt.geom_signal(sig, REGION, n_bins=500, agg="mean")
    expected = G.signal_steps(sig, REGION, n_bins=500, agg="mean")
    pd.testing.assert_frame_equal(layer.data[["start", "end", "value"]], expected)


def test_signal_fill_and_outline():
    sig = pygt.example_signal()
    assert isinstance(pgt.geom_signal(sig, REGION), geom_rect)
    outline = pgt.geom_signal(sig, REGION, fill=False)
    assert isinstance(outline, geom_path)
    # The example signal has one gap, so the outline has two runs
    assert outline.data["_path"].nunique() == 2


def test_signal_overlay_bins_each_sample():
    sig = pygt.example_signal()
    two = pd.concat([sig.assign(sample="a"), sig.assign(sample="b", signal=sig.signal * 2)])
    layer = pgt.geom_signal(two, REGION, mapping=aes(fill="sample"), n_bins=300)
    assert set(layer.data["sample"]) == {"a", "b"}
    a = layer.data.query("sample == 'a'")["value"].to_numpy()
    b = layer.data.query("sample == 'b'")["value"].to_numpy()
    np.testing.assert_allclose(b, 2 * a)
    assert "fill" not in layer.aes_params  # mapped, so not fixed


def test_signal_by_keeps_facet_column():
    sig = pygt.example_signal()
    two = pd.concat([sig.assign(cond="wt"), sig.assign(cond="ko")])
    p = ggplot() + pgt.geom_signal(two, REGION, by=["cond"]) + facet_wrap("cond")
    p += pgt.scale_x_genome(REGION)
    fig = p.draw()
    assert len(fig.axes) == 2


def test_signal_data_empty_region():
    sig = pygt.example_signal()
    out = pgt.signal_data(sig, "chr9:1-100")
    assert list(out.columns) == ["start", "end", "value"] and out.empty


# Intervals


def test_intervals_pack_and_label():
    peaks = pygt.example_peaks()
    frame = pgt.interval_data(peaks, REGION, label_by="name")
    assert {"xmin", "xmax", "row", "y", "ymin", "ymax"} <= set(frame.columns)
    # boxes in one row never overlap
    for _, row in frame.groupby("row"):
        row = row.sort_values("xmin")
        assert (row["xmin"].to_numpy()[1:] >= row["xmax"].to_numpy()[:-1]).all()
    layers = pgt.geom_intervals(peaks, REGION, aes(fill="cre"), label_by="name")
    assert any(isinstance(lyr, geom_text) for lyr in layers)
    assert "fill" not in layers[0].aes_params


def test_intervals_min_width_keeps_motifs_visible():
    motifs = pygt.example_motif_hits()
    frame = pgt.interval_data(motifs, REGION, min_width_px=3, panel_width=5, dpi=100)
    assert ((frame["xmax"] - frame["xmin"]) >= 3 * REGION.length / 500 - 1e-6).all()


def test_intervals_strand_arrowheads():
    motifs = pygt.example_motif_hits()
    big = motifs.assign(end=motifs["start"] + 3000)
    layers = pgt.geom_intervals(big, REGION, strand="strand")
    shapes = {
        lyr.aes_params["shape"]
        for lyr in layers
        if "shape" in getattr(lyr, "aes_params", {})
    }
    assert shapes <= {">", "<"} and shapes


# Genes


def test_genes_from_dataframe():
    spans, parts = pgt.gene_data(pygt.example_genes(), REGION)
    assert set(parts["transcript"]) <= set(spans["transcript"])
    full = parts[parts["full"]]
    utr = parts[~parts["full"]]
    assert ((full["ymax"] - full["ymin"]) > (utr["ymax"] - utr["ymin"]).max()).all()
    assert (spans["start"] >= REGION.start).all() and (spans["end"] <= REGION.end).all()


def test_genes_longest():
    spans, _ = pgt.gene_data(pygt.example_genes(), REGION, transcripts="longest")
    assert spans["name"].is_unique


def test_genes_empty_region_still_gives_axis():
    layers = pgt.geom_genes(pygt.example_genes(), "chr9:1-1000")
    assert render(ggplot() + layers + pgt.scale_x_genome("chr9:1-1000")) > 0


# Layout


def tracks():
    xs = pgt.scale_x_genome(REGION)
    sig = ggplot() + pgt.geom_signal(pygt.example_signal(), REGION) + xs + pgt.theme_track()
    cre = ggplot() + pgt.geom_intervals(pygt.example_peaks(), REGION, label_by="name") + xs
    gen = ggplot() + pgt.geom_genes(pygt.example_genes(), REGION) + xs + labs(y="Genes")
    return sig + labs(y="Msn2 ChEC\n(rpm)"), cre + pgt.theme_track(), gen + pgt.theme_track()


def test_stack_tracks_aligns_panels():
    comp = pgt.stack_tracks(*tracks())
    fig = comp.draw()
    FigureCanvasAgg(fig).draw()
    xs = {(round(a.get_window_extent().x0, 3), round(a.get_window_extent().x1, 3)) for a in fig.axes}
    assert len(xs) == 1
    # only the bottom track keeps x tick labels
    shown = [any(t.get_visible() and t.get_text() for t in a.get_xticklabels()) for a in fig.axes]
    assert shown == [False, False, True]


def test_stack_tracks_heights():
    try:
        from plotnine.composition import plot_layout  # noqa: F401
    except ImportError:
        with pytest.raises(PlotnineError, match="plot_layout"):
            pgt.stack_tracks(*tracks(), heights=[2, 1, 1])
        return
    comp = pgt.stack_tracks(*tracks(), heights=[2, 1, 1])
    buf = io.BytesIO()
    comp.save(buf, format="png", width=6, height=5, dpi=72, verbose=False)
    fig = comp.draw()
    FigureCanvasAgg(fig).draw()
    h = [a.get_window_extent().height for a in fig.axes]
    assert h[0] / h[1] == pytest.approx(2, rel=0.02)


def test_stack_tracks_arguments():
    with pytest.raises(ValueError):
        pgt.stack_tracks(tracks()[0])
    with pytest.raises(ValueError):
        pgt.stack_tracks(*tracks(), heights=[1, 1])


@pytest.mark.parametrize("fmt", ["png", "pdf", "svg"])
def test_formats(fmt):
    comp = pgt.stack_tracks(*tracks())
    buf = io.BytesIO()
    comp.save(buf, format=fmt, width=6, height=4, verbose=False)
    assert buf.getbuffer().nbytes > 1000
