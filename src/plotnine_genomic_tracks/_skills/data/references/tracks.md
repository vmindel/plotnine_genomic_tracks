# Tracks

Every track function returns plotnine layers; add them to `ggplot()`
with `+`, together with `pgt.scale_x_genome(region)`.

## pgt.scale_x_genome

```python
pgt.scale_x_genome(region, unit=None, name=True, expand=False, ylim=None)
```

- Returns a scale **and** `coord_cartesian(xlim=region)`. Adding another
  coord (`coord_cartesian`, `coord_flip`) afterwards replaces the zoom.
- Ticks in bp/kb/Mb chosen from the region length; `unit="kb"` forces
  one. `name=True` titles the axis `"chr2 (kb)"`.
- Zooms rather than filters: data outside the region is kept, so
  `annotate("rect", ...)` highlights and `geom_vline` need no clipping.
- `ylim=(0, 60)`: y-range of a **signal** plot, cut at the panel edge
  (bins above it are drawn cut, not dropped). Never on an interval or gene
  plot: it overrides their row layout and squashes the features.

## pgt.geom_signal

```python
pgt.geom_signal(df, region, value="signal", mapping=None, *, by=None,
                cols=("chrom", "start", "end"), n_bins=2000, agg="max",
                fill=True, color="#3b6ea8", size=0.4, **kwargs)
```

- Draws each row as a step over its real extent. Gaps between rows and
  NaN values stay empty; they are never interpolated.
- More rows than `n_bins`, or overlapping rows: aggregated onto `n_bins`
  equal bins across the region. `agg="max"` (default) keeps narrow peaks;
  `"mean"` is overlap-weighted; `"min"`. `n_bins=None` never aggregates
  and raises on overlapping rows.
- `fill=True`: filled area, colour from `color`. `fill=False`: outline,
  colour from `color`, width `size`. `**kwargs` go to `geom_rect` or
  `geom_path` (e.g. `alpha=0.5`).
- **Overlays**: `mapping=aes(color="sample")` with `fill=False`, or
  `aes(fill="sample")` with `alpha`. Every column named in `mapping` is
  binned separately and the y-scale is shared. Map only categorical
  columns: a numeric column makes one series per distinct value.
- **Facets**: pass the facet column in `by=["condition"]`, then
  `facet_wrap("condition", ncol=1)`. Without `by` the column is not in
  the layer data and faceting fails.
- No y-range label in IGV style: the y-axis is a normal axis.

## pgt.geom_intervals

```python
pgt.geom_intervals(df, region, mapping=None, *, cols=..., label_by=None,
                   strand=None, pack=True, height=0.6, min_width_px=1.0,
                   fill="#595959", color=None, fontsize=7,
                   panel_width=6.0, dpi=100, **kwargs)
```

- Any BED-like table: peaks, CREs, motif hits. Overlapping features are
  packed into rows; with `label_by`, rows leave room for the labels.
- Colour is a normal aesthetic: `aes(fill="cre")` with
  `scale_fill_manual({...})` (dict keeps colours stable across figures),
  or `aes(fill="score")` with `scale_fill_cmap("viridis")`. Legends work.
- `strand="strand"` draws white arrowheads on features at least 6 px
  wide; narrower ones get none.
- `min_width_px` widens short features around their centre; raise it
  (e.g. 4-6) for motif hits in a wide window.
- The result includes the track's y scale (rows, no ticks). One interval
  or gene track per plot.
- All columns of `df` are kept in the layer, so faceting on any of them
  works.

## pgt.geom_genes

```python
pgt.geom_genes(genes, region, *, transcripts="all", cols=...,
               feature="feature", transcript="transcript_id",
               name="gene_name", strand="strand", color="#333333",
               height=0.6, utr_ratio=0.5, show_names=True, pack=True,
               fontsize=7, panel_width=6.0, dpi=100)
```

- `genes`: a UCSC assembly name (`"hg38"`, `"mm10"`, `"sacCer3"`), a GTF
  path, or a DataFrame with one row per exon/CDS (`pygt.read_gtf`
  output). Genome names download from UCSC on first use and are cached;
  on offline cluster nodes run `pygt download mm10` on a login node
  first.
- Intron line, exons, CDS at full `height`, UTRs at `height * utr_ratio`,
  strand arrowheads, names under each transcript.
- `transcripts="longest"`: one transcript per gene, usually what a figure
  needs.
- Not facetable: the layer data does not keep extra columns.

## Resolution assumptions

plotnine fixes layer data before it knows the panel size. So:

- signal is binned to `n_bins=2000` regardless of panel width;
- packing, `min_width_px` and arrowhead spacing assume a panel
  `panel_width=6` inches wide at `dpi=100`.

For a narrower figure where labels collide, pass the real panel width
(roughly figure width minus 1.5 inches for labels and legend).
