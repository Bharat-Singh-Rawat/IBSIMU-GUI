"""Shared matplotlib utilities for tkinter embedding."""
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image


def make_canvas(parent, figsize=(5, 3.3), dpi=100):
    """Create a Figure + FigureCanvasTkAgg pair."""
    fig = Figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=parent)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    return fig, ax, canvas


def placeholder(ax, text="No data", title=""):
    ax.clear()
    ax.text(0.5, 0.5, text, ha="center", va="center",
            transform=ax.transAxes, color="gray")
    if title:
        ax.set_title(title)


def draw_rms_ellipse(ax, eps, alpha, beta):
    """Draw Courant-Snyder RMS ellipse on an axes."""
    if eps <= 0 or beta <= 0:
        return
    t = np.linspace(0, 2 * np.pi, 200)
    sb = np.sqrt(beta)
    y = np.sqrt(eps) * sb * np.cos(t)
    yp = np.sqrt(eps) * (-alpha / sb * np.cos(t) + 1 / sb * np.sin(t))
    ax.plot(y, yp, "r-", lw=1.5, alpha=0.8, label="RMS ellipse")


def capture_composite(figs, layout="1+2"):
    """Capture multiple figures into a composite PIL image.
    layout='1+2': first fig on top, remaining side-by-side below."""
    imgs = []
    for fig in figs:
        fig.canvas.draw()
        buf = fig.canvas.buffer_rgba()
        w, h = fig.canvas.get_width_height()
        imgs.append(Image.frombuffer("RGBA", (w, h), buf).copy())
    if not imgs:
        return None
    top = imgs[0]
    tw = top.width
    if len(imgs) >= 3:
        bh = max(imgs[1].height, imgs[2].height)
        hw = tw // 2
        comp = Image.new("RGBA", (tw, top.height + bh), (255, 255, 255, 255))
        comp.paste(top, (0, 0))
        comp.paste(imgs[1].resize((hw, bh), Image.LANCZOS), (0, top.height))
        comp.paste(imgs[2].resize((tw - hw, bh), Image.LANCZOS), (hw, top.height))
    elif len(imgs) == 2:
        comp = Image.new("RGBA", (tw, top.height + imgs[1].height), (255, 255, 255, 255))
        comp.paste(top, (0, 0))
        comp.paste(imgs[1].resize((tw, imgs[1].height), Image.LANCZOS), (0, top.height))
    else:
        comp = top
    return comp.convert("RGB")
