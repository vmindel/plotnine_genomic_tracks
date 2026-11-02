# Stacking tracks

`stack_tracks(*plots, heights=None)` stacks track plots top to bottom.
Panels line up exactly whatever the width of each track's y-axis labels,
and only the bottom track keeps its x-axis.

```python
pgt.stack_tracks(signal, cres, motifs, genes, heights=[1.4, 0.8, 0.6, 0.8])
```

Save the result with `.save()`:

```python
pgt.stack_tracks(signal, genes).save("locus.pdf", width=7, height=4)
```

## Track heights need plotnine's development version

`heights` uses plotnine's `plot_layout`, which is in plotnine's
development branch and not in 0.15.8. Without it, `stack_tracks` stacks
the tracks at equal heights, and asking for `heights` raises an error
that says so. To get heights now:

```bash
pip install git+https://github.com/has2k1/plotnine
```

## Saving

Use the composition's `.save()`. Saving a figure that was already drawn
again through matplotlib (`fig = comp.draw(); fig.savefig(...)` more than
once) loses the track heights in current plotnine; `.save()` lays the
figure out afresh each time and is not affected.

## Labels and legends

`theme_track()` writes track labels (the y-axis title, set with
`labs(y=...)`) horizontally beside the panel and keeps legends compact so
they fit beside short tracks. Each track keeps its own legend.
