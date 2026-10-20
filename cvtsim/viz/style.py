"""Visual system: monochrome, one red accent, dot-matrix display type (Doto) over a
monospaced body (Space Mono). Both fonts are SIL OFL, vendored from google/fonts."""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

FONTS = pathlib.Path(__file__).with_name("fonts")

THEMES = {
    "dark": dict(BLACK="#000000", PANEL="#0d0d0d", GRID="#262626", MUTED="#7a7a7a", WHITE="#f2f2f2",
                 RED="#d71921", SERIES=["#f2f2f2", "#d71921", "#9a9a9a", "#555555"], CMAP=(0.3, 0.95)),
    "light": dict(BLACK="#f4f4f1", PANEL="#ffffff", GRID="#d6d6d2", MUTED="#6b6b6b", WHITE="#111111",
                  RED="#c4141c", SERIES=["#111111", "#c4141c", "#7a7a7a", "#b0b0b0"], CMAP=(0.0, 0.75)),
}
# Active palette, rebound by use(). WHITE is the foreground ink in either theme.
BLACK = PANEL = GRID = MUTED = WHITE = RED = ""
SERIES = []
CMAP_RANGE = (0.0, 1.0)
theme = None


def use(name=None):
    """Register fonts and apply a theme ('dark' or 'light'); no argument keeps the current one."""
    global BLACK, PANEL, GRID, MUTED, WHITE, RED, SERIES, CMAP_RANGE, theme
    name = name or theme or "dark"
    if name == theme:
        return
    if theme is None:
        for f in FONTS.glob("*.ttf"):
            font_manager.fontManager.addfont(str(f))
    theme = name
    t = THEMES[name]
    BLACK, PANEL, GRID, MUTED, WHITE, RED = (t[k] for k in ("BLACK", "PANEL", "GRID", "MUTED", "WHITE", "RED"))
    SERIES, CMAP_RANGE = list(t["SERIES"]), t["CMAP"]
    plt.rcParams.update({
        "figure.facecolor": BLACK, "axes.facecolor": PANEL, "savefig.facecolor": BLACK,
        "axes.edgecolor": GRID, "axes.labelcolor": MUTED, "axes.titlecolor": WHITE,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": (0, (1, 3)),
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": WHITE,
        "font.family": ["Space Mono", "DejaVu Sans Mono"], "font.size": 8.5,
        "axes.titlesize": 11, "axes.titlelocation": "left", "axes.titlepad": 10,
        "axes.labelsize": 8, "legend.fontsize": 7.5, "legend.frameon": False,
        "axes.prop_cycle": plt.cycler(color=SERIES),
        "lines.linewidth": 1.3, "lines.solid_capstyle": "round",
        "figure.dpi": 110, "savefig.dpi": 160,
    })


def cmap():
    """Perceptually uniform sequential map (magma segment) that stays visible on the theme."""
    import numpy as np
    from matplotlib.colors import ListedColormap
    return ListedColormap(plt.get_cmap("magma")(np.linspace(*CMAP_RANGE, 256)), "cvtsim")


def gradient_line(ax, x, y, c, lw=2.0, vmin=None, vmax=None):
    """Polyline coloured by per-point value c through the perceptual map. Returns the collection."""
    import numpy as np
    from matplotlib.collections import LineCollection
    pts = np.column_stack([x, y])[:, None, :]
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
    c = np.asarray(c, dtype=float)
    lc = LineCollection(seg, cmap=cmap(), linewidths=lw, capstyle="round")
    lc.set_array(0.5 * (c[:-1] + c[1:]))
    lc.set_clim(c.min() if vmin is None else vmin, c.max() if vmax is None else vmax)
    ax.add_collection(lc)
    ax.autoscale_view()
    return lc


def colorbar(fig, mappable, ax, label=""):
    cb = fig.colorbar(mappable, ax=ax, pad=0.02, fraction=0.05)
    cb.outline.set_visible(False)
    cb.set_label(label, color=MUTED, fontsize=8)
    cb.ax.tick_params(colors=MUTED, labelsize=7.5)
    return cb


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
