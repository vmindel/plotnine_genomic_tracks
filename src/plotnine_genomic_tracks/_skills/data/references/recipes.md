# Recipes

## From pymisha

`gextract` returns `chrom`, `start`, `end`, one column per expression
(named by `colnames`) and `intervalID`; pass the column as `value`.
`gscreen` and interval sets return `chrom`, `start`, `end`.

```python
import pymisha as pm
import pygenomic_tracks as pygt
import plotnine_genomic_tracks as pgt
from plotnine import aes, ggplot, labs

region = pygt.gene_region("mm10", "Gata2", flank=50_000)
scope = pm.gintervals(region.chrom, region.start, region.end)
chec = pm.gextract("chec.msn2", intervals=scope, iterator=50, colnames=["msn2"])
peaks = pm.gscreen("chec.msn2 > 20", intervals=scope)

x = pgt.scale_x_genome(region)
fig = pgt.stack_tracks(
    ggplot() + pgt.geom_signal(chec, region, value="msn2")
    + x + pgt.theme_track() + labs(y="Msn2"),
    ggplot() + pgt.geom_intervals(peaks, region) + x + pgt.theme_track()
    + labs(y="Peaks"),
    ggplot() + pgt.geom_genes("mm10", region, transcripts="longest")
    + x + pgt.theme_track() + labs(y="Genes"),
)
fig.save("gata2.pdf", width=7, height=4)
```

`pygt.gene_region(genome, gene, flank=...)` returns a `pygt.Region`
around a gene.

## Samples overlaid in one track

```python
import pandas as pd
from plotnine import scale_color_manual
both = pd.concat([msn2.assign(sample="Msn2"), msn4.assign(sample="Msn4")])
p = (ggplot()
     + pgt.geom_signal(both, region, mapping=aes(color="sample"), fill=False)
     + scale_color_manual({"Msn2": "#3b6ea8", "Msn4": "#c0504d"})
     + pgt.scale_x_genome(region) + pgt.theme_track() + labs(y="ChEC"))
```

## Conditions as facets, one shared scale

```python
from plotnine import facet_wrap
both = pd.concat([wt.assign(cond="WT"), ko.assign(cond="KO")])
p = (ggplot()
     + pgt.geom_signal(both, region, by=["cond"], agg="mean")
     + facet_wrap("cond", ncol=1)
     + pgt.scale_x_genome(region) + pgt.theme_track())
```

## CREs coloured by class, motif hits by score

```python
from plotnine import scale_fill_cmap, scale_fill_manual
cres_p = (ggplot()
          + pgt.geom_intervals(cres, region, aes(fill="cre"), label_by="name")
          + scale_fill_manual({"promoter": "#c0504d", "enhancer": "#3b6ea8"})
          + x + pgt.theme_track() + labs(y="CREs", fill=""))
motif_p = (ggplot()
           + pgt.geom_intervals(motifs, region, aes(fill="score"),
                                strand="strand", min_width_px=4)
           + scale_fill_cmap("Greys") + x + pgt.theme_track())
```

## Dense signal in a vector figure

A 1 Mb window at 10 bp resolution makes a heavy PDF. With
plotnine-extensions, rasterise only the signal layer:

```python
from plotnine_extensions import rasterise
p = ggplot() + rasterise(pgt.geom_signal(chec, region), dpi=300) + x
```

## Example data

`pygt.example_region()`, `pygt.example_signal()`, `pygt.example_peaks()`,
`pygt.example_motif_hits()` and `pygt.example_genes()` give a small
synthetic locus for trying things out without files.
