"""
Genome-browser tracks as plotnine layers

Built on pygenomic_tracks: the same binning, packing and gene models as
its matplotlib tracks, drawn as plotnine layers that facet, compose and
theme like any other.
"""

from .axis import scale_x_genome
from .features import gene_data, geom_genes, geom_intervals, interval_data
from .layout import stack_tracks, theme_track
from .signal import geom_signal, signal_data

__version__ = "0.1.0"

__all__ = [
    "gene_data",
    "geom_genes",
    "geom_intervals",
    "geom_signal",
    "interval_data",
    "scale_x_genome",
    "signal_data",
    "stack_tracks",
    "theme_track",
]
