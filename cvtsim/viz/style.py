"""Visual system: monochrome, one red accent, dot-matrix display type (Doto) over a
monospaced body (Space Mono). Both fonts are SIL OFL, vendored from google/fonts."""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

FONTS = pathlib.Path(__file__).with_name("fonts")

BLACK = "#000000"
PANEL = "#0d0d0d"
GRID = "#262626"
MUTED = "#7a7a7a"
WHITE = "#f2f2f2"
RED = "#d71921"
# series order: white, red, then greys (never more than four series per panel)
SERIES = [WHITE, RED, "#9a9a9a", "#555555"]

_ready = False


def use():
    """Register fonts and apply the theme (idempotent)."""
    global _ready
    if _ready:
        return
    for f in FONTS.glob("*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({
        "figure.facecolor": BLACK, "axes.facecolor": PANEL, "savefig.facecolor": BLACK,
        "axes.edgecolor": GRID, "axes.labelcolor": MUTED, "axes.titlecolor": WHITE,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": (0, (1, 3)),
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": WHITE,
        "font.family": "Space Mono", "font.size": 8.5,
        "axes.titlesize": 11, "axes.titlelocation": "left", "axes.titlepad": 10,
        "axes.labelsize": 8, "legend.fontsize": 7.5, "legend.frameon": False,
        "axes.prop_cycle": plt.cycler(color=SERIES),
        "lines.linewidth": 1.3, "lines.solid_capstyle": "round",
        "figure.dpi": 110, "savefig.dpi": 160,
    })
    _ready = True


def title(ax, text):
    """Panel title in the dot-matrix face, upper case."""
    ax.set_title(text.upper(), fontfamily="Doto", fontweight=900, fontsize=12)


def header(fig, text, sub=""):
    fig.text(0.012, 0.985, text.upper(), fontfamily="Doto", fontweight=900,
             fontsize=22, color=WHITE, va="top")
    if sub:
        fig.text(0.012, 0.945, sub, fontsize=8, color=MUTED, va="top")


def card(ax, label, value, unit="", accent=False):
    """Metric tile: small grey label, large dot-matrix value."""
    ax.set_axis_off()
    ax.set_facecolor(PANEL)
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=PANEL, zorder=-1))
    ax.text(0.06, 0.78, label.upper(), fontsize=7, color=MUTED, transform=ax.transAxes)
    ax.text(0.06, 0.18, value, fontfamily="Doto", fontweight=900, fontsize=24,
            color=RED if accent else WHITE, transform=ax.transAxes)
    if unit:
        ax.text(0.94, 0.2, unit, fontsize=7.5, color=MUTED, ha="right", transform=ax.transAxes)


def ema(y, alpha):
    """Exponential moving average, as used for smoothed training curves."""
    import numpy as np
    out = np.empty_like(y, dtype=float)
    acc = y[0]
    for i, v in enumerate(y):
        acc = alpha * acc + (1 - alpha) * v
        out[i] = acc
    return out
