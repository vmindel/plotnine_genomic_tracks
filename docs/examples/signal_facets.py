import pandas as pd
import pygenomic_tracks as pygt
from plotnine import facet_wrap, ggplot, labs

import plotnine_genomic_tracks as pgt

region = "chr2:1,215,000-1,285,000"
chec = pygt.example_signal()
conditions = pd.concat(
    [chec.assign(condition="wild type"), chec.assign(condition="Δmsn2", signal=chec["signal"] * 0.3)]
)

# `by` bins each condition separately and keeps the column to facet on
p = (
    ggplot()
    + pgt.geom_signal(conditions, region, by=["condition"], agg="mean")
    + facet_wrap("condition", ncol=1)
    + pgt.scale_x_genome(region)
    + pgt.theme_track()
    + labs(y="ChEC")
)
