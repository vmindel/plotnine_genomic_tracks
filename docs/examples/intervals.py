import pygenomic_tracks as pygt
from plotnine import aes, ggplot, labs, scale_fill_manual

import plotnine_genomic_tracks as pgt

region = pygt.example_region()

p = (
    ggplot()
    + pgt.geom_intervals(pygt.example_peaks(), region, aes(fill="cre"), label_by="name")
    + scale_fill_manual({"promoter": "#c0504d", "enhancer": "#3b6ea8", "CTCF": "#8064a2"})
    + pgt.scale_x_genome(region)
    + pgt.theme_track()
    + labs(y="CREs", fill="")
)
