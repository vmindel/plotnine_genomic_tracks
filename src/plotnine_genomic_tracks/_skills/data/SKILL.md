---
name: plotnine-genomic-tracks
description: Use for genome-browser figures built with plotnine (ggplot-style, layers added with +) via plotnine_genomic_tracks (import plotnine_genomic_tracks as pgt). Triggers on plotting a genomic region or locus in plotnine; signal / coverage / ChEC-seq / ChIP-seq / ATAC tracks, peaks, CREs, motif hits or gene models as plotnine layers; stacking tracks on one genomic x-axis; faceting tracks by sample or condition; Gviz or plotgardener-style figures in plotnine. For the matplotlib API (pygt.signal(df, ax, region)) use the pygenomic-tracks skill instead. Router skill: read the relevant file under references/ before writing code.
---

# plotnine_genomic_tracks

Genome-browser tracks as plotnine layers, built on pygenomic_tracks: the
same binning, row packing and gene models as its matplotlib tracks, drawn
as layers that facet, take `scale_fill_*`/`scale_color_*` and legends,
and stack on one genomic x-axis. This file only routes; the reference
files hold signatures, defaults and the mistakes that give a wrong figure
without an error. Read the file for the task before writing code.

## Task to reference file

| task | read |
|---|---|
| **one track**: `pgt.geom_signal`, `pgt.geom_intervals`, `pgt.geom_genes`, `pgt.scale_x_genome`; arguments, defaults, colour, facets, overlays | `references/tracks.md` |
| **a multi-track figure**: `pgt.stack_tracks`, heights, `pgt.theme_track`, shared y-range, highlights, several loci, saving | `references/layout.md` |
| **worked examples**: pymisha `gextract`/`gscreen` to tracks, sample overlays, conditions as facets, motif zoom, rasterising dense signal | `references/recipes.md` |
| "**why is my track empty / wrong?**", or before handing over a figure | `references/pitfalls.md` |

## Conventions

- `import plotnine_genomic_tracks as pgt` and, for example data and
  readers, `import pygenomic_tracks as pygt`.
- Track functions take **the table and the region**, not plot-level data:
  `ggplot() + pgt.geom_signal(df, region, value="msn2")`. They ignore
  `ggplot(df, aes(...))`.
- `region` is `"chr2:1,200,000-1,350,000"`, `("chr2", 1200000, 1350000)`
  or `pygt.Region`; 0-based, half-open. Pass the **same** region to every
  track function and to `pgt.scale_x_genome(region)` in every plot.
- Tables need `chrom`, `start`, `end` (0-based, half-open) plus value
  columns; rename with `cols=("chr", "from", "to")`. pymisha's `gextract`
  and `gscreen` output goes in unchanged. `chr2` and `2` match each other.
- One track per plot; stack plots with
  `pgt.stack_tracks(p1, p2, ..., heights=[...])` and save with `.save()`.
- `pgt.theme_track()` gives the genome-browser look; set the track label
  with `labs(y="...")`.
- To see what a track will draw, call `pgt.signal_data`,
  `pgt.interval_data` or `pgt.gene_data` with the same arguments.
