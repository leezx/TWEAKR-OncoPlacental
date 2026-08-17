#!/usr/bin/env python3
"""
Step 9: render the 3 Fetal-Placenta-Adult ternary maps (Track A/B/C) as
static PNGs. Standard barycentric-to-Cartesian projection: Adult=(0,0),
Placenta=(1,0), Fetal=(0.5, sqrt(3)/2); x = Placenta + 0.5*Fetal,
y = (sqrt(3)/2)*Fetal.

Usage: python3 plot_ternary.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

OUT_DIR = "/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental/results/09_developmental_ternary"
SQRT3_2 = np.sqrt(3) / 2

TRACKS = [
    ("track_A_gut_specific_coords.tsv", "Track A -- Gut-specific\n(Fetal colon vs Placenta trophoblast vs Adult colon)", "ternary_track_A_gut_specific.png"),
    ("track_B_pantissue_coords.tsv", "Track B -- Pan-tissue\n(Fetal-somatic [7 HDMA organs] vs Placenta trophoblast vs Adult [GTEx+HPA])", "ternary_track_B_pantissue.png"),
    ("track_C_hcl_coords.tsv", "Track C -- HCL independent validation (GSE134355)\n(Fetal intestine vs Placenta [n=1] vs Adult intestine, single atlas)", "ternary_track_C_hcl.png"),
]


def to_xy(fetal, placenta, adult):
    x = placenta + 0.5 * fetal
    y = SQRT3_2 * fetal
    return x, y


def draw_triangle(ax):
    verts = [(0, 0), (1, 0), (0.5, SQRT3_2)]
    ax.add_patch(Polygon(verts, closed=True, fill=False, edgecolor="black", linewidth=1.2, zorder=5))
    ax.text(0.5, SQRT3_2 + 0.03, "Fetal", ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax.text(-0.06, -0.04, "Adult", ha="right", va="top", fontsize=13, fontweight="bold")
    ax.text(1.06, -0.04, "Placenta", ha="left", va="top", fontsize=13, fontweight="bold")
    # gridlines at 25/50/75% for each axis (light guide)
    for frac in (0.25, 0.5, 0.75):
        # lines of constant Fetal
        x0, y0 = to_xy(frac, 1 - frac, 0)
        x1, y1 = to_xy(frac, 0, 1 - frac)
        ax.plot([x0, x1], [y0, y1], color="gray", lw=0.4, alpha=0.5, zorder=1)
        # lines of constant Placenta
        x0, y0 = to_xy(1 - frac, frac, 0)
        x1, y1 = to_xy(0, frac, 1 - frac)
        ax.plot([x0, x1], [y0, y1], color="gray", lw=0.4, alpha=0.5, zorder=1)
        # lines of constant Adult
        x0, y0 = to_xy(1 - frac, 0, frac)
        x1, y1 = to_xy(0, 1 - frac, frac)
        ax.plot([x0, x1], [y0, y1], color="gray", lw=0.4, alpha=0.5, zorder=1)


def plot_track(path, title, out_name):
    # keep_default_na=False: the "marker" column's empty-string rows (not
    # a known marker gene) must round-trip as "" not NaN -- pandas'
    # default NA-sniffing on read_csv otherwise turns every empty field
    # into NaN, and `df["marker"] != ""` then evaluates True for EVERY
    # row (NaN != "" is always True), mislabeling the whole gene universe
    # as markers. Real bug, caught by looking at the first rendered plot.
    df = pd.read_csv(f"{OUT_DIR}/{path}", sep="\t", index_col=0, keep_default_na=False)
    x, y = to_xy(df["Fetal"].astype(float).values, df["Placenta"].astype(float).values, df["Adult"].astype(float).values)

    fig, ax = plt.subplots(figsize=(8, 7.2))
    draw_triangle(ax)
    ax.hexbin(x, y, gridsize=60, cmap="Blues", mincnt=1, zorder=2, alpha=0.85)

    markers = df[df["marker"] != ""] if "marker" in df.columns else df.iloc[0:0]
    if len(markers):
        mx, my = to_xy(markers["Fetal"].values, markers["Placenta"].values, markers["Adult"].values)
        ax.scatter(mx, my, color="red", s=36, zorder=10, edgecolor="black", linewidth=0.6)
        # Marker points frequently collapse near-exactly onto each other
        # (e.g. several placental hormone genes all near the Placenta
        # vertex) -- stack their labels vertically instead of letting
        # xytext overlap into unreadable jumbled text (real readability
        # bug, caught by looking at the first rendered Track B/C plots).
        order = np.argsort(-my)  # stable top-to-bottom stacking
        placed = []
        for rank, i in enumerate(order):
            gx, gy, gname = mx[i], my[i], markers["marker"].values[i]
            y_off = 6 + 13 * sum(1 for px, py in placed if abs(px - gx) < 0.03 and abs(py - gy) < 0.03)
            placed.append((gx, gy))
            ax.annotate(gname, (gx, gy), textcoords="offset points", xytext=(6, y_off),
                        fontsize=9, fontweight="bold", color="darkred", zorder=11)

    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.1, SQRT3_2 + 0.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{title}\n(n={len(df)} genes)", fontsize=12)
    fig.tight_layout()
    out_path = f"{OUT_DIR}/{out_name}"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    for path, title, out_name in TRACKS:
        plot_track(path, title, out_name)


if __name__ == "__main__":
    main()
