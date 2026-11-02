import pygenomic_tracks as pygt
from plotnine import aes, ggplot, labs, scale_fill_cmap

import plotnine_genomic_tracks as pgt

region = "chr2:1,262,000-1,276,000"  # zoom in on E1 and E2

p = (
    ggplot()
    + pgt.geom_intervals(
        pygt.example_motif_hits(),
        region,
        aes(fill="score"),
        strand="strand",
        label_by="motif",
        min_width_px=6,
    )
    + scale_fill_cmap("viridis")
    + pgt.scale_x_genome(region)
    + pgt.theme_track()
    + labs(y="Motifs")
)
