import pandas as pd
import pygenomic_tracks as pygt
from plotnine import aes, ggplot, labs, scale_color_manual

import plotnine_genomic_tracks as pgt

region = "chr2:1,215,000-1,285,000"
chec = pygt.example_signal()
samples = pd.concat(
    [chec.assign(sample="Msn2"), chec.assign(sample="Msn4", signal=chec["signal"] * 0.6)]
)

# Each sample is binned separately; the y scale is shared
p = (
    ggplot()
    + pgt.geom_signal(samples, region, mapping=aes(color="sample"), fill=False, size=0.5)
    + scale_color_manual(["#3b6ea8", "#c0504d"])
    + pgt.scale_x_genome(region)
    + pgt.theme_track()
    + labs(y="ChEC", color="")
)
