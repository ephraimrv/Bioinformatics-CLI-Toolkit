"""
Figure 3 Generator: Promoter Region Architecture Schematic

Draws a high-resolution, vector-based structural comparison of two promoter
regions to illustrate the absence of regulatory motifs in foreign loci.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example usage:
    $ python3 generate_figure3.py
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.3"

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- FIGURE PARAMETERS & MAGIC NUMBERS ---
TRACK_Y_ANCESTRAL = 2.0
TRACK_Y_FOREIGN = -1.0
MOTIF_START = -85
MOTIF_LENGTH = 15

COLOR_MOTIF = "#ce93d8"  # Orchid purple — Regulatory Element
COLOR_ATG_ANCESTRAL = "#4caf50"  # Green
COLOR_ATG_FOREIGN = "#ef5350"  # Red
COLOR_DNA_BACKBONE = "#555555"
COLOR_RULER = "#888888"


def main():
    fig, ax = plt.subplots(figsize=(10, 5))

    # Clean canvas
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

    # Add Panel Identifier (Top Left)
    ax.text(-170, 3.0, "B", fontsize=16, weight="bold", fontfamily="sans-serif")

    # ---------------------------------------------------------
    # TRACK 1: Ancestral Lactobin Promoter
    # ---------------------------------------------------------
    ax.plot(
        [-150, 0],
        [TRACK_Y_ANCESTRAL, TRACK_Y_ANCESTRAL],
        color=COLOR_DNA_BACKBONE,
        lw=2,
    )

    # Directional Start Codon Box (Ancestral)
    arrow_atg1 = patches.FancyArrow(
        0,
        TRACK_Y_ANCESTRAL,
        20,
        0,
        width=0.3,
        head_width=0.5,
        head_length=5,
        color=COLOR_ATG_ANCESTRAL,
        alpha=0.9,
        length_includes_head=True,
    )
    ax.add_patch(arrow_atg1)
    ax.text(
        5, TRACK_Y_ANCESTRAL + 0.3, "ctg1_68\n(Lactobin A)", fontsize=10, weight="bold"
    )
    ax.text(3, TRACK_Y_ANCESTRAL - 0.35, "ATG", fontsize=8, family="monospace")

    rect_motif = patches.Rectangle(
        (MOTIF_START, TRACK_Y_ANCESTRAL - 0.1),
        MOTIF_LENGTH,
        0.2,
        color=COLOR_MOTIF,
        zorder=3,
    )
    ax.add_patch(rect_motif)
    ax.text(
        MOTIF_START - 10,
        TRACK_Y_ANCESTRAL + 0.2,
        "Motif 1 Operator",
        fontsize=9,
        color="#7b1fa2",
        weight="bold",
    )

    # Scale Ruler
    ruler_y1 = TRACK_Y_ANCESTRAL - 0.7
    ax.plot([-150, 0], [ruler_y1, ruler_y1], color=COLOR_RULER, lw=1, linestyle="--")
    for tick in [-150, -100, -50, 0]:
        ax.plot(
            [tick, tick], [ruler_y1 - 0.05, ruler_y1 + 0.05], color=COLOR_RULER, lw=1
        )
        ax.text(tick - 5, ruler_y1 - 0.2, f"{tick} bp", fontsize=8, color="#666666")

    # ---------------------------------------------------------
    # TRACK 2: Foreign Blp Promoter
    # ---------------------------------------------------------
    ax.plot(
        [-150, 0], [TRACK_Y_FOREIGN, TRACK_Y_FOREIGN], color=COLOR_DNA_BACKBONE, lw=2
    )

    # Directional Start Codon Box (Foreign)
    arrow_atg2 = patches.FancyArrow(
        0,
        TRACK_Y_FOREIGN,
        20,
        0,
        width=0.3,
        head_width=0.5,
        head_length=5,
        color=COLOR_ATG_FOREIGN,
        alpha=0.9,
        length_includes_head=True,
    )
    ax.add_patch(arrow_atg2)
    ax.text(
        5, TRACK_Y_FOREIGN + 0.3, "ctg1_50\n(Blp locus)", fontsize=10, weight="bold"
    )
    ax.text(3, TRACK_Y_FOREIGN - 0.35, "ATG", fontsize=8, family="monospace")

    ax.text(
        -120,
        TRACK_Y_FOREIGN + 0.2,
        "[ No Regulatory Motifs Discovered ]",
        fontsize=10,
        color="#d32f2f",
        style="italic",
    )
    ax.text(
        -120,
        TRACK_Y_FOREIGN - 0.3,
        "Unregulated, highly A/T-rich promoter space",
        fontsize=9,
        color="#666666",
    )

    # Scale Ruler
    ruler_y2 = TRACK_Y_FOREIGN - 0.7
    ax.plot([-150, 0], [ruler_y2, ruler_y2], color=COLOR_RULER, lw=1, linestyle="--")
    for tick in [-150, -100, -50, 0]:
        ax.plot(
            [tick, tick], [ruler_y2 - 0.05, ruler_y2 + 0.05], color=COLOR_RULER, lw=1
        )
        ax.text(tick - 5, ruler_y2 - 0.2, f"{tick} bp", fontsize=8, color="#666666")

    # Set Field Boundaries
    ax.set_xlim(-170, 50)
    ax.set_ylim(-2.5, 3.5)

    # ---------------------------------------------------------
    # SAFE EXPORT
    # ---------------------------------------------------------
    pdf_out = "Figure3_PanelB_Promoters.pdf"
    tiff_out = "Figure3_PanelB_Promoters.tiff"

    try:
        plt.savefig(pdf_out, bbox_inches="tight")
        plt.savefig(tiff_out, format="tiff", dpi=300, bbox_inches="tight")
        print(f"[✓] Panel B schematic successfully generated!")
        print(f"    -> {pdf_out} (Vector)")
        print(f"    -> {tiff_out} (300 DPI Raster)")
    except OSError as e:
        sys.exit(f"\n[!] Error saving schematic: {e}")
    finally:
        plt.close(fig)


if __name__ == "__main__":
    main()
