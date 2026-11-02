from __future__ import annotations

from typing import TYPE_CHECKING

from plotnine import coord_cartesian, scale_x_continuous
from pygenomic_tracks import as_region
from pygenomic_tracks.geometry import genomic_unit

if TYPE_CHECKING:
    from pygenomic_tracks._coords import RegionLike

__all__ = ["scale_x_genome"]


def scale_x_genome(
    region: RegionLike,
    unit: str | None = None,
    name: str | bool = True,
    expand: bool = False,
) -> list:
    """
    Genomic x-axis for a region

    Zooms the panel to the region and labels positions in bp, kb or Mb.
    Data outside the region is kept, so layers such as `geom_vline`
    highlights need no clipping.

    Parameters
    ----------
    region : Region, str or tuple
        ``"chr2:1,200,000-1,350,000"``, ``("chr2", 1200000, 1350000)`` or a
        `pygenomic_tracks.Region`; 0-based, half-open.
    unit : {"bp", "kb", "Mb"}, optional
        Unit of the labels. By default the one that suits the region
        length, as in pygenomic_tracks.
    name : str or bool
        Axis title. `True` gives e.g. ``"chr2 (kb)"``, `False` none.
    expand : bool
        Pad the panel beyond the region (and the y data range). Genome
        browsers do not; the default is `False`.

    Returns
    -------
    list
        A scale and a coordinate system; add both with ``+``.
    """
    region = as_region(region)
    units = {"bp": 1.0, "kb": 1e3, "Mb": 1e6}
    if unit is None:
        unit, factor = genomic_unit(region.length)
    elif unit in units:
        factor = units[unit]
    else:
        raise ValueError(f"unit must be one of {list(units)}, got {unit!r}")

    def labels(breaks):
        return [f"{b / factor:,.3f}".rstrip("0").rstrip(".") for b in breaks]

    if name is True:
        title = f"{region.chrom} ({unit})"
    elif name is False:
        title = None
    else:
        title = name

    return [
        scale_x_continuous(labels=labels, name=title),
        coord_cartesian(xlim=(region.start, region.end), expand=expand),
    ]
