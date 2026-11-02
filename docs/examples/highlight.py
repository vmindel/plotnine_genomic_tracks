import numpy as np
import pygenomic_tracks as pygt
from plotnine import annotate, ggplot, labs

import plotnine_genomic_tracks as pgt

region = pygt.example_region()
x = pgt.scale_x_genome(region)

# Ordinary plotnine layers mix in: here a highlight added to every track
enhancers = annotate(
    "rect", xmin=1_260_000, xmax=1_274_000, ymin=-np.inf, ymax=np.inf,
    fill="#f2c14e", alpha=0.25,
)

signal = (
    ggplot() + enhancers + pgt.geom_signal(pygt.example_signal(), region)
    + x + pgt.theme_track() + labs(y="ChEC-seq")
)
genes = (
    ggplot() + enhancers + pgt.geom_genes(pygt.example_genes(), region, transcripts="longest")
    + x + pgt.theme_track() + labs(y="Genes")
)

p = pgt.stack_tracks(signal, genes, heights=[1.5, 1])
