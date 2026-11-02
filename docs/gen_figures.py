"""
Run every docs example and save the figure it defines

Each file in docs/examples defines `p`; the docs include the same files
verbatim, so every figure comes from the code shown next to it.

    python docs/gen_figures.py
"""

import runpy
import warnings
from pathlib import Path

HERE = Path(__file__).parent
SIZES = {
    "quickstart": (7, 5),
    "highlight": (7, 3.6),
    "signal_facets": (7, 3.6),
    "genes": (7, 2.4),
    "intervals": (7, 2.2),
    "motifs": (7, 2.6),
    "signal_overlay": (7, 2.4),
}


def main():
    (HERE / "img").mkdir(exist_ok=True)
    for path in sorted((HERE / "examples").glob("*.py")):
        width, height = SIZES.get(path.stem, (7, 3))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = runpy.run_path(str(path))["p"]
            p.save(
                HERE / "img" / f"{path.stem}.png",
                width=width, height=height, dpi=130, verbose=False,
            )
        print(f"docs/img/{path.stem}.png")


if __name__ == "__main__":
    main()
