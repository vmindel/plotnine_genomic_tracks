# Figure layout

## pgt.stack_tracks

```python
fig = pgt.stack_tracks(signal, cres, genes, heights=[1.5, 0.8, 1.0])
fig.save("locus.pdf", width=7, height=4)
```

- Plots top to bottom. Panels line up exactly whatever each track's
  y-label width; only the bottom plot keeps its x-axis.
- Every plot must carry the same `pgt.scale_x_genome(region)`; stacking
  aligns panels, not data ranges.
- `heights` needs plotnine's `plot_layout`, which is only in plotnine's
  development branch (`pip install git+https://github.com/has2k1/plotnine`),
  not in 0.15.8. On 0.15.8 omit `heights` (equal heights) or install the
  development version; with `heights` it raises a `PlotnineError` that says
  so. Check with `from plotnine.composition import plot_layout`.

## Saving

- Call `.save(path, width=, height=)` on the composition (or plot)
  directly; width and height are inches.
- Do not save through matplotlib (`fig = comp.draw(); fig.savefig(...)`):
  track heights drift from what `heights` asked for. On plotnine 0.15.8,
  `comp.draw()` followed by `comp.save()` on a side-by-side (`|`)
  composition raises `KeyError`. Build the composition and save it once.

## Look

- `pgt.theme_track(base_size=9)`: classic theme, horizontal track label
  beside the panel (`labs(y="Msn2 ChEC")`), compact legends. Add it to
  each plot. Feature tracks hide their y ticks through their scale, so the
  order of `+ theme_track()` does not matter.
- Each track keeps its own legend at its right.

## Comparable signal tracks

Separate plots have separate y-scales: two samples in two tracks are not
comparable by eye unless they share a range:

```python
x = pgt.scale_x_genome(region)
x_sig = pgt.scale_x_genome(region, ylim=(0, 60))   # signal plots only
```

Use `x_sig` on every signal plot and `x` on interval and gene plots.
Overlays inside one plot share a scale automatically.

## Highlights and annotations

Plots are ordinary plotnine plots; add the same layer to every track:

```python
import numpy as np
from plotnine import annotate
hl = annotate("rect", xmin=1_261_000, xmax=1_274_000, ymin=-np.inf,
              ymax=np.inf, fill="#f2c14e", alpha=0.25)
signal = ggplot() + hl + pgt.geom_signal(chec, region) + x + pgt.theme_track()
```

Add the highlight before the track so it is drawn underneath.

## Several loci

One figure per locus is simplest. Side by side, compose stacks with `|`:

```python
def locus(region):
    x = pgt.scale_x_genome(region)
    x_sig = pgt.scale_x_genome(region, ylim=(0, 60))
    return pgt.stack_tracks(
        ggplot() + pgt.geom_signal(chec, region) + x_sig + pgt.theme_track(),
        ggplot() + pgt.geom_genes("mm10", region, transcripts="longest")
        + x + pgt.theme_track(),
    )

(locus(r1) | locus(r2)).save("loci.pdf", width=10, height=3)
```

Give signal tracks the same `ylim` across loci so peaks compare.
