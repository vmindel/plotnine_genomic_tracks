# Pitfalls

Mistakes that give a wrong or empty figure without an error. Check these
before handing over a figure.

## Empty track

- **Region and data on different chromosomes or assemblies.** `chr2` and
  `2` match each other; anything else (`chrM` vs `MT`, hg19 data in an
  hg38 region) gives an empty track. Check with
  `pgt.signal_data(df, region, value=...)` or `pgt.interval_data(df, region)`:
  an empty table means nothing overlaps.
- **Plot-level data.** `ggplot(df) + pgt.geom_signal(...)` ignores `df`;
  the table must be the track function's first argument.
- **Wrong value column.** `geom_signal` reads `value="signal"` by default;
  pymisha names it after `colnames`. A missing column raises `KeyError`,
  but an all-NaN column draws nothing.
- **1-based coordinates** (GTF, R, IGV) are off by one base; invisible at
  track scale but wrong in tables. `pygt.read_gtf` converts GTFs.

## Wrong-looking track

- **Different regions in one figure.** The region in a track function
  selects and clips data; `scale_x_genome(region)` sets the axis. If they
  differ, data is cut where the axis is not, or the panel shows empty
  space. Use one `region` variable everywhere.
- **`ylim` on an interval or gene plot** squashes the rows; use it only
  on signal plots (see `references/layout.md`).
- **Adding a coord after `scale_x_genome`** (`coord_cartesian(ylim=...)`,
  `coord_flip()`) replaces the zoom; the panel shows the whole data range.
  Use `scale_x_genome(region, ylim=...)` instead.
- **Signal tracks compared by eye with independent y-scales.** Each plot
  scales to its own maximum; share `ylim` across tracks and loci.
- **A numeric column in `mapping` of `geom_signal`** splits the signal into
  one series per distinct value. Map categorical columns only.
- **Two interval or gene tracks in one plot**: each sets the y scale; the
  last one wins and the other is drawn in the wrong rows. plotnine only
  prints "Scale for 'y' is already present". One feature track per plot.
- **Colliding labels in a narrow figure**: packing assumes a 6 inch
  panel; pass `panel_width=` near the real panel width.
- **Labels near the panel edge are clipped**, as in a genome browser.

## Figure assembly

- **`heights` on plotnine 0.15.8** raises `PlotnineError`: `plot_layout`
  is only in plotnine's development branch.
- **`fig = comp.draw()` then `fig.savefig(...)`**: track heights drift
  from what `heights` asked for. On 0.15.8, `comp.draw()` followed by
  `comp.save()` on a side-by-side (`|`) composition raises `KeyError`.
  Build the composition and call `.save()` once.
- **Faceting a signal track without `by=`** raises "At least one layer
  must contain all variables used for facetting"; pass `by=["column"]`.
- **Category colours that change between figures**: with
  `scale_fill_manual`, pass a dict `{category: colour}`, not a list.
