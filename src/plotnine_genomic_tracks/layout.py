from __future__ import annotations

from typing import TYPE_CHECKING

from plotnine import element_blank, element_line, element_text, theme, theme_classic
from plotnine.composition import Stack
from plotnine.exceptions import PlotnineError

if TYPE_CHECKING:
    from collections.abc import Sequence

    from plotnine import ggplot

__all__ = ["stack_tracks", "theme_track"]


def theme_track(base_size: float = 9, base_family: str | None = None) -> theme:
    """
    Genome-browser look: no grid, track label beside the panel

    Parameters
    ----------
    base_size : float
        Base font size in points.
    base_family : str, optional
        Font family.
    """
    return theme_classic(base_size=base_size, base_family=base_family) + theme(
        # Horizontal track labels. Left-aligned: plotnine reserves too
        # little room for right-aligned horizontal y titles.
        axis_title_y=element_text(rotation=0, ha="left", va="center"),
        axis_line=element_line(size=0.5),
        # Compact legends, so a legend fits beside a short track
        legend_key_size=base_size - 1,
        legend_text=element_text(size=base_size - 2),
        legend_title=element_text(size=base_size - 1),
    )


_HIDE_X = theme(
    axis_text_x=element_blank(),
    axis_title_x=element_blank(),
    axis_ticks_x=element_blank(),
    axis_line_x=element_blank(),
)


def stack_tracks(*plots: ggplot, heights: Sequence[float] | None = None):
    """
    Stack track plots on one genomic x-axis

    Panels line up exactly whatever their y-axis labels; only the bottom
    track keeps its x-axis.

    Parameters
    ----------
    *plots : ggplot
        Track plots, top to bottom. Give each the same `scale_x_genome`.
    heights : sequence of float, optional
        Relative track heights. Needs plotnine's ``plot_layout``, which is
        in plotnine's development branch and not in 0.15.8; without it,
        tracks share the height equally.

    Returns
    -------
    plotnine composition
        Draw, show or save it like a plot.
    """
    if len(plots) < 2:
        raise ValueError("stack_tracks needs at least two plots")
    items = [p + _HIDE_X for p in plots[:-1]] + [plots[-1]]
    stack = Stack(items)
    if heights is None:
        return stack
    if len(heights) != len(plots):
        raise ValueError(
            f"got {len(heights)} heights for {len(plots)} tracks"
        )
    try:
        from plotnine.composition import plot_layout
    except ImportError as err:
        raise PlotnineError(
            "Track heights need plotnine's plot_layout, which this plotnine "
            "version does not have. Install the development version "
            "(pip install git+https://github.com/has2k1/plotnine) or stack "
            "the tracks without heights."
        ) from err
    return stack + plot_layout(heights=list(heights))
