import pygenomic_tracks as pygt
from plotnine import aes, ggplot, labs, scale_fill_cmap

import plotnine_genomic_tracks as pgt

region = pygt.example_region()  # chr2:1,200,000-1,350,000
x = pgt.scale_x_genome(region)

signal = (
    ggplot()
    + pgt.geom_signal(pygt.example_signal(), region)
    + x + pgt.theme_track() + labs(y="ChEC-seq")
)
cres = (
    ggplot()
    + pgt.geom_intervals(pygt.example_peaks(), region, aes(fill="cre"), label_by="name")
    + x + pgt.theme_track() + labs(y="CREs")
)
motifs = (
    ggplot()
    + pgt.geom_intervals(pygt.example_motif_hits(), region, aes(fill="score"))
    + scale_fill_cmap("Greys")
    + x + pgt.theme_track() + labs(y="Motifs")
)
genes = (
    ggplot()
    + pgt.geom_genes(pygt.example_genes(), region, transcripts="longest")
    + x + pgt.theme_track() + labs(y="Genes")
)

p = pgt.stack_tracks(signal, cres, motifs, genes, heights=[1.4, 0.8, 0.6, 0.8])
