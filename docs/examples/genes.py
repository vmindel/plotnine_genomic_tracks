import pygenomic_tracks as pygt
from plotnine import ggplot, labs

import plotnine_genomic_tracks as pgt

region = pygt.example_region()

# A DataFrame here; "mm10", "hg38" or a GTF path work the same way
p = (
    ggplot()
    + pgt.geom_genes(pygt.example_genes(), region)
    + pgt.scale_x_genome(region)
    + pgt.theme_track()
    + labs(y="Genes")
)
