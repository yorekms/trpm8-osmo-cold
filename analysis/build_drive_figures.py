#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy>=1.26",
#   "matplotlib>=3.8",
#   "pillow>=10.0",
# ]
# ///
"""
Compose 100% original, publication-grade, copyright-clean Figures 1 and 2 for
bioRxiv / journal submission (ZERO third-party BioRender assets):

  - analysis/outputs/fig1_corneal_nerve_interface.png
      Panel A: Pure cropped mouse corneal subbasal nerve whorl microscopy (Yorek stock)
      Panel B: Original vector anatomical & osmotic-barrier comparison
               (Avascular Corneal Epithelium vs. Vascularized Cutaneous Epidermis)
      Panel C: Intraepithelial free nerve terminals triple-stain cross-section (30 um scale; Yorek stock)

  - analysis/outputs/fig2_mechanistic_circuit.png
      100% original vector biophysical schematic of the coupled electro-metabolic
      TRPM8 - Nav1.8 - Cav1.2/1.3 - Mitochondrial MICU1/MCU/TMBIM5 H+ oscillator:
      Panel A: Basal Poised State (resting Vm ~ -45 mV, MICU1 gate closed)
      Panel B: Osmo-Thermal Burst & L-Type Amplification ([Ca2+]i -> 515 nM)
      Panel C: Mitochondrial H+ Negative-Feedback Reset & Clinical/Genetic Lesions
"""

import os
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse, Polygon, PathPatch
from matplotlib.path import Path as MplPath
from PIL import Image, ImageOps


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LEGACY_LINUX_DRIVE_DIR = Path("/home/yorekms/.gemini/antigravity/brain/598ce7af-ddbd-4f27-9c36-189ac5f71b07/scratch/drive/figures")


def resolve_microscopy_asset(canonical_name: str, drive_Fallbacks: list[Path]) -> Path:
    """Locate a source microscopy image across env override, local assets/, synced Google Drive, or Linux scratch."""
    candidates: list[Path] = []
    env_dir = os.environ.get("TRPM8_ASSETS_DIR")
    if env_dir:
        candidates.append(Path(env_dir) / canonical_name)
    candidates.append(SCRIPT_DIR / "assets" / canonical_name)
    candidates.extend(drive_Fallbacks)
    candidates.append(LEGACY_LINUX_DRIVE_DIR / canonical_name)
    for p in candidates:
        if p.is_file():
            return p
    searched = "\n  - ".join(str(p) for p in candidates)
    raise FileNotFoundError(f"Could not find microscopy asset '{canonical_name}'. Searched:\n  - {searched}")


def trim_whitespace(img: Image.Image, pad: int = 10) -> Image.Image:
    rgb = img.convert("RGB")
    bg = rgb.getpixel((2, 2))
    if min(bg) > 240:
        inv = ImageOps.invert(rgb)
        bbox = inv.getbbox()
        if bbox:
            x0, y0, x1, y1 = bbox
            w, h = rgb.size
            return rgb.crop((max(0, x0 - pad), max(0, y0 - pad), min(w, x1 + pad), min(h, y1 + pad)))
    return rgb


def draw_cornea_vs_skin_panel(ax) -> None:
    """Original vector anatomical & water-activity comparison: Avascular Cornea vs. Vascularized Skin."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # ---- LEFT HALF: AVASCULAR CORNEA (x = 0.5..43.5) ----
    ax.add_patch(FancyBboxPatch((0.5, 3.0), 43.0, 93.0, boxstyle="round,pad=0.6", facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.2))
    ax.text(22.0, 92.6, "Avascular Cornea (Benchmark)", ha="center", va="center",
            fontsize=8.6, fontweight="bold", color="#0369a1")

    # Evaporating Tear Film
    ax.add_patch(mpatches.Rectangle((2.5, 72.5), 39.0, 12.0, facecolor="#bae6fd", edgecolor="#38bdf8", lw=1.0, alpha=0.88))
    ax.text(22.0, 81.0, "Evaporating Tear Film", ha="center", va="center", fontsize=8.0, fontweight="bold", color="#0c4a6e")
    ax.text(22.0, 75.8, "300 → 450+ mOsm (>340 excites; Parra 2014)", ha="center", va="center", fontsize=6.4, fontweight="bold", color="#0369a1")
    for ex in [8.5, 22.0, 35.5]:
        ax.annotate("", xy=(ex, 89.4), xytext=(ex, 84.7),
                    arrowprops=dict(arrowstyle="->", color="#0284c7", lw=1.5))

    # Corneal Epithelium (thin, non-keratinized, 5 cell layers)
    ax.add_patch(mpatches.Rectangle((2.5, 35.0), 39.0, 37.5, facecolor="#e0f2fe", edgecolor="#7dd3fc", lw=0.9))
    for cy in [42.5, 50.5, 58.5, 66.5]:
        ax.plot([3.0, 41.0], [cy, cy], color="#93c5fd", lw=0.65, alpha=0.65)
        for cx in np.linspace(6.8, 37.2, 5):
            ax.add_patch(Ellipse((cx, cy - 3.2), 3.8, 2.5, facecolor="#60a5fa", edgecolor="none", alpha=0.40))
    ax.text(22.0, 39.8, "Non-Keratinized Epithelium (~40 µm)", ha="center", va="center",
            fontsize=7.1, fontweight="bold", color="#1e3a8a",
            bbox=dict(boxstyle="round,pad=0.16", facecolor="#e0f2fe", edgecolor="none", alpha=0.93), zorder=8)

    # Subbasal Plexus + Avascular Stroma
    ax.add_patch(mpatches.Rectangle((2.5, 5.5), 39.0, 29.5, facecolor="#f1f5f9", edgecolor="#cbd5e1", lw=0.9))
    ax.text(22.0, 22.5, "Avascular Stroma\n(No Blood Vessels | ΔFlow = 0)", ha="center", va="center",
            fontsize=7.6, fontweight="bold", color="#334155")
    ax.text(22.0, 10.8, "Direct aqueous osmo-thermal access\nZero vascular washout confound", ha="center", va="center",
            fontsize=7.0, style="italic", color="#0f766e")

    # Subbasal nerve trunk + vertical piercing intraepithelial TRPM8+ terminals
    ax.plot([4.5, 39.5], [34.8, 34.8], color="#dc2626", lw=2.4, zorder=5)
    for nx in [10.0, 19.0, 28.0, 35.0]:
        xs = [nx, nx - 1.0, nx + 0.8, nx - 0.5]
        ys = [34.8, 48.0, 61.0, 71.8]
        ax.plot(xs, ys, color="#dc2626", lw=1.9, zorder=6)
        ax.scatter([xs[-1]], [ys[-1]], s=36, color="#9333ea", edgecolor="white", lw=0.65, zorder=7)
    ax.text(22.0, 63.0, "TRPM8⁺ Complex Bulbar Endings\n(Mito-Packed; Brock 2001, 2018)", ha="center", va="center", fontsize=6.5, fontweight="bold",
            color="#991b1b", bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#fecaca", alpha=0.95), zorder=8)

    # ---- RIGHT HALF: VASCULARIZED SKIN (x = 46.0..99.5) ----
    ax.add_patch(FancyBboxPatch((46.0, 3.0), 53.5, 93.0, boxstyle="round,pad=0.6", facecolor="#fffbeb", edgecolor="#fcd34d", lw=1.2))
    ax.text(72.75, 92.6, "Vascularized Skin (SC Gradient & SG2 Frontier)", ha="center", va="center",
            fontsize=8.2, fontweight="bold", color="#92400e")

    # Dry Air / High VPD
    ax.add_patch(mpatches.Rectangle((48.0, 79.0), 49.5, 6.2, facecolor="#fef3c7", edgecolor="#fbbf24", lw=0.8))
    ax.text(72.75, 82.1, r"Dry Air (20–50% RH) → TEWL | $a_w \approx 0.20$–$0.50$", ha="center", va="center",
            fontsize=6.9, fontweight="bold", color="#92400e")
    for ex in [55, 72.75, 90.5]:
        ax.annotate("", xy=(ex, 89.4), xytext=(ex, 85.4),
                    arrowprops=dict(arrowstyle="->", color="#d97706", lw=1.5))

    # Upper / Mid Stratum Corneum (NMF Zone)
    ax.add_patch(mpatches.Rectangle((48.0, 68.0), 49.5, 11.0, facecolor="#fde68a", edgecolor="#d97706", lw=0.9))
    for sy in [70.8, 73.5, 76.2]:
        ax.plot([48.5, 97.0], [sy, sy], color="#b45309", lw=0.55, alpha=0.55)
    ax.text(72.75, 73.5, r"Mid SC (NMF Zone): $a_w\ 0.50\rightarrow 0.95$ (~3–12 Osm/kg)", ha="center", va="center",
            fontsize=6.7, fontweight="bold", color="#78350f",
            bbox=dict(boxstyle="round,pad=0.14", facecolor="#fffbeb", edgecolor="none", alpha=0.92))

    # Basal Cornified Layer (Stratum Compactum)
    ax.add_patch(mpatches.Rectangle((48.0, 58.6), 49.5, 9.4, facecolor="#fed7aa", edgecolor="#c2410c", lw=0.9))
    for sy in [61.7, 64.9]:
        ax.plot([48.5, 97.0], [sy, sy], color="#9a3412", lw=0.5, alpha=0.5)
    ax.text(72.75, 63.3, r"Basal SC (Compactum): low NMF | inner $a_w \approx 0.995$", ha="center", va="center",
            fontsize=6.6, fontweight="bold", color="#7c2d12",
            bbox=dict(boxstyle="round,pad=0.14", facecolor="#fff7ed", edgecolor="none", alpha=0.93))

    # Living Epidermis: Stratum Granulosum (upper band) + Spinosum/Basale (lower band)
    ax.add_patch(mpatches.Rectangle((48.0, 46.8), 49.5, 11.8, facecolor="#ffedd5", edgecolor="#fb923c", lw=0.8))
    ax.add_patch(mpatches.Rectangle((48.0, 31.8), 49.5, 15.0, facecolor="#fff7ed", edgecolor="#fdba74", lw=0.8))

    # SG2 Tight-Junction Seal (dashed emerald line + junction beads)
    ax.plot([48.5, 97.0], [56.8, 56.8], color="#0f766e", lw=1.8, ls="--", zorder=6)
    for tx in np.linspace(51.5, 94.0, 9):
        ax.add_patch(mpatches.Rectangle((tx - 0.7, 56.1), 1.4, 1.4, facecolor="#0f766e", edgecolor="white", lw=0.4, zorder=7))
    ax.text(72.75, 52.6, "SG2 Tight Junctions (Prune Nerves) • NFAT5/TAUT", ha="center", va="center",
            fontsize=6.6, fontweight="bold", color="#0f766e",
            bbox=dict(boxstyle="round,pad=0.14", facecolor="#f0fdfa", edgecolor="#99f6e4", alpha=0.96), zorder=8)

    # Cutaneous TRPM8+ C-fiber terminals ascending to JUST BELOW SG2 tight junctions
    ax.plot([51, 94.5], [31.8, 31.8], color="#7c3aed", lw=2.2, zorder=5)
    for nx in [57.5, 72.75, 88.0]:
        xs = [nx, nx - 1.1, nx + 0.9, nx - 0.3]
        ys = [31.8, 41.0, 49.5, 55.5]
        ax.plot(xs, ys, color="#7c3aed", lw=1.8, zorder=6)
        ax.scatter([xs[-1]], [ys[-1]], s=28, color="#9333ea", edgecolor="white", lw=0.6, zorder=7)

    # Osmotic leverage callout right in living epidermis
    ax.text(72.75, 40.4, "Cutaneous TRPM8⁺ Endings (Dhaka 2008; Takahashi 2019)\n" + r"99.8% SC shield $\rightarrow \Delta a_w = -0.001 \rightarrow +56$ mOsm (+89 winter)", ha="center", va="center",
            fontsize=6.4, fontweight="bold", color="#5b21b6",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#ddd6fe", alpha=0.96), zorder=9)

    # Perfused Dermis with Capillary Loop & Cav1.2 Arteriole
    ax.add_patch(mpatches.Rectangle((48.0, 5.5), 49.5, 26.3, facecolor="#fee2e2", edgecolor="#fca5a5", lw=0.9))
    theta = np.linspace(0, np.pi, 40)
    ax.plot(63.0 + 5.2 * np.cos(theta), 21.5 + 8.0 * np.sin(theta), color="#dc2626", lw=3.5, alpha=0.82, zorder=3)
    ax.plot(82.5 + 5.2 * np.cos(theta), 21.5 + 8.0 * np.sin(theta), color="#ef4444", lw=3.5, alpha=0.82, zorder=3)
    ax.add_patch(FancyBboxPatch((51, 17.2), 43.5, 5.8, boxstyle="round,pad=0.25,rounding_size=1.8",
                                facecolor="#b91c1c", edgecolor="#7f1d1d", lw=0.9, zorder=4))
    ax.text(72.75, 20.1, "Caᵥ1.2 Dermal Arteriole (~37 °C)", ha="center", va="center",
            fontsize=7.4, fontweight="bold", color="white", zorder=5)
    ax.text(72.75, 10.6, "Clamps ~290 mOsm/kg + Warm Blood\n(DHP vasodilation confounds neural block; §3.4)", ha="center", va="center",
            fontsize=6.8, fontweight="bold", color="#991b1b")


def build_figure_1() -> Path:
    """Build 100% copyright-free Figure 1 (pure Yorek stock microscopy + original vector diagram)."""
    cornea_path = resolve_microscopy_asset(
        "mouse_cornea.jpeg",
        [
            SCRIPT_DIR.parent.parent / "figures" / "mouse cornea.jpeg",
            SCRIPT_DIR.parent.parent / "figures" / "mouse_cornea.jpeg",
        ],
    )
    triple_path = resolve_microscopy_asset(
        "cornea_triple.jpg",
        [
            SCRIPT_DIR.parent.parent / "figures" / "cornea_triple.jpg",
            SCRIPT_DIR.parent.parent.parent.parent / "Photos 2016" / "Cornea figure - triple.jpg",
        ],
    )
    raw_cornea = Image.open(cornea_path).convert("RGB")
    img_whorl = raw_cornea.crop((2078, 2460, 4742, 4112))
    img_triple = trim_whitespace(Image.open(triple_path), pad=6)

    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(10.8, 8.6), dpi=240, facecolor="white")
    gs = gridspec.GridSpec(
        2, 2,
        height_ratios=[1.22, 0.78],
        width_ratios=[0.74, 1.38],
        hspace=0.20,
        wspace=0.05,
        left=0.025, right=0.975, top=0.905, bottom=0.03,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    ax_a.imshow(img_whorl)
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    for s in ax_a.spines.values():
        s.set_edgecolor("#94a3b8")
        s.set_linewidth(1.0)
    ax_a.set_title("A   Corneal Subbasal Nerve Whorl (C / Aδ)",
                   loc="left", fontsize=9.8, fontweight="bold", color="#0f172a", pad=7)

    draw_cornea_vs_skin_panel(ax_b)
    ax_b.set_title("B   Anatomical Osmotic-Barrier Contrast: Cornea vs. Skin",
                   loc="left", fontsize=9.8, fontweight="bold", color="#0f172a", pad=7)

    ax_c.imshow(img_triple)
    ax_c.set_xticks([])
    ax_c.set_yticks([])
    for s in ax_c.spines.values():
        s.set_edgecolor("#94a3b8")
        s.set_linewidth(1.0)
    ax_c.set_title("C   Intraepithelial Corneal Free Nerve Terminals Piercing to Tear Surface (Fluorescein / TO-PRO-3 / βIII-Tubulin; 30 µm Scale)",
                   loc="left", fontsize=9.8, fontweight="bold", color="#0f172a", pad=7)

    fig.suptitle(
        "Figure 1. Anatomical Interface for Peripheral TRPM8 Osmosensation: Corneal Benchmark vs. Cutaneous Barrier",
        fontsize=11.8,
        fontweight="bold",
        color="#0f172a",
        y=0.968,
    )

    out_path = OUT_DIR / "fig1_corneal_nerve_interface.png"
    fig.savefig(out_path, dpi=240, facecolor="white")
    plt.close(fig)
    return out_path


# ==================== FIGURE 2: FULL-PAGE PRINT-READABLE 3-ROW VECTOR CIRCUIT ====================

def draw_channel_wide(ax, cx: float, cy: float, color: str, edge: str, open_pore: bool, label: str, sublabel: str = ""):
    """Draw a large, crisp print-readable two-lobe transmembrane ion channel on a 3:1 aspect panel."""
    gap = 1.65 if open_pore else 0.28
    w = 2.15
    h = 18.5
    ax.add_patch(FancyBboxPatch((cx - gap / 2 - w, cy - h / 2), w, h, boxstyle="round,pad=0.35,rounding_size=1.1",
                                facecolor=color, edgecolor=edge, lw=1.6, zorder=5))
    ax.add_patch(FancyBboxPatch((cx + gap / 2, cy - h / 2), w, h, boxstyle="round,pad=0.35,rounding_size=1.1",
                                facecolor=color, edgecolor=edge, lw=1.6, zorder=5))
    if open_pore:
        ax.add_patch(mpatches.Rectangle((cx - gap / 2, cy - h / 2 - 0.6), gap, h + 1.2,
                                        facecolor="white", edgecolor="none", alpha=0.90, zorder=4))
    ax.text(cx, cy + h / 2 + 2.4, label, ha="center", va="bottom", fontsize=10.0, fontweight="bold", color=edge, zorder=8)
    if sublabel:
        ax.text(cx, cy - h / 2 - 3.0, sublabel, ha="center", va="top", fontsize=8.6, fontweight="bold", color="#334155", zorder=8)


def draw_mitochondrion_wide(ax, cx: float, cy: float, rx: float, ry: float, ca_loaded: bool, gate_open: bool, h_release: bool):
    """Draw a large, print-readable mitochondrion with cristae, MICU1/MCU pore, and TMBIM5."""
    # Outer membrane
    ax.add_patch(Ellipse((cx, cy), 2 * rx, 2 * ry, facecolor="#fed7aa", edgecolor="#9a3412", lw=1.8, zorder=3))
    # Inner matrix
    ax.add_patch(Ellipse((cx, cy), 2 * (rx - 1.5), 2 * (ry - 2.8), facecolor="#fff7ed", edgecolor="#c2410c", lw=1.2, zorder=4))

    # Crisp cristae strictly inside inner ellipse
    ri_x = rx - 2.0
    ri_y = ry - 3.6
    xs_c = np.linspace(cx - ri_x + 3.2, cx + ri_x - 3.2, 7)
    for idx, xk in enumerate(xs_c):
        dx_norm = (xk - cx) / ri_x
        y_edge = ri_y * np.sqrt(max(0.0, 1.0 - dx_norm * dx_norm)) - 0.8
        if idx % 2 == 0:
            ax.plot([xk, xk], [cy - y_edge, cy + 2.0], color="#c2410c", lw=1.3, alpha=0.68, zorder=5)
        else:
            ax.plot([xk, xk], [cy + y_edge, cy - 2.0], color="#c2410c", lw=1.3, alpha=0.68, zorder=5)

    # MCU + MICU1 complex on upper-left rim
    gx, gy = cx - 7.6, cy + ry - 1.5
    gate_col = "#15803d" if gate_open else "#64748b"
    ax.add_patch(FancyBboxPatch((gx - 3.4, gy - 4.0), 6.8, 8.0, boxstyle="round,pad=0.35",
                                facecolor=gate_col, edgecolor="white", lw=1.1, zorder=7))
    ax.text(gx, gy, "MCU", ha="center", va="center", fontsize=8.6, fontweight="bold", color="white", zorder=8)
    ax.add_patch(FancyBboxPatch((gx - 6.8, gy + 5.2), 13.6, 7.2, boxstyle="round,pad=0.3",
                                facecolor="#dcfce7" if gate_open else "#f1f5f9", edgecolor=gate_col, lw=1.2, zorder=7))
    ax.text(gx, gy + 8.8, "MICU1 " + ("OPEN" if gate_open else "CLOSED"), ha="center", va="center",
            fontsize=8.1, fontweight="bold", color=gate_col, zorder=8)

    # TMBIM5 Ca2+/H+ exchanger on upper-right rim
    tx, ty = cx + 7.8, cy + ry - 1.5
    ex_col = "#dc2626" if h_release else "#64748b"
    ax.add_patch(FancyBboxPatch((tx - 5.2, ty - 4.0), 10.4, 8.0, boxstyle="round,pad=0.35",
                                facecolor=ex_col, edgecolor="white", lw=1.1, zorder=7))
    ax.text(tx, ty, "TMBIM5", ha="center", va="center", fontsize=8.4, fontweight="bold", color="white", zorder=8)

    # Matrix status label
    if ca_loaded:
        ax.text(cx, cy - 2.2, "Matrix Ca²⁺ ↑  (ΔΨm)", ha="center", va="center", fontsize=9.0, fontweight="bold",
                color="#9a3412", bbox=dict(boxstyle="round,pad=0.22", facecolor="#ffedd5", edgecolor="#fb923c", alpha=0.96), zorder=9)
    else:
        ax.text(cx, cy - 2.2, "Mitochondrion (ΔΨm intact)", ha="center", va="center", fontsize=8.8, fontweight="bold",
                color="#7c2d12", bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="none", alpha=0.90), zorder=9)


def draw_wide_stage_frame(ax, bg_ext: str, stage_title: str, badge_col: str, ext_caption: str, readout_text: str):
    """Draw a full-width horizontal nerve terminal panel sized for crisp print readability."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # Outer extracellular frame
    ax.add_patch(FancyBboxPatch((0.6, 1.5), 98.8, 97.0, boxstyle="round,pad=0.4", facecolor=bg_ext, edgecolor="#cbd5e1", lw=1.3, zorder=0))

    # Cytosolic nerve terminal pocket spanning full width below membrane (y = 5..66)
    ax.add_patch(FancyBboxPatch((2.2, 5.0), 95.6, 61.0, boxstyle="round,pad=0.4,rounding_size=3.0",
                                facecolor="white", edgecolor="#475569", lw=3.0, zorder=1))
    ax.add_patch(FancyBboxPatch((3.0, 6.4), 94.0, 58.2, boxstyle="round,pad=0.3,rounding_size=2.4",
                                facecolor="#fcfdff", edgecolor="#cbd5e1", lw=1.0, zorder=2))

    # Left-aligned stage header pill cleanly separated from channel labels below
    ax.text(2.8, 92.8, stage_title, ha="left", va="center", fontsize=9.0, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.32", facecolor=badge_col, edgecolor="none"), zorder=10)
    # Right-aligned extracellular condition text
    ax.text(97.4, 92.8, ext_caption, ha="right", va="center", fontsize=8.5, fontweight="bold", color="#1e293b", zorder=10)

    # Crisp bottom physiology readout bar inside terminal
    ax.add_patch(FancyBboxPatch((4.2, 7.8), 91.6, 11.8, boxstyle="round,pad=0.3",
                                facecolor="#f1f5f9", edgecolor="#94a3b8", lw=0.9, zorder=6))
    ax.text(50.0, 13.7, readout_text, ha="center", va="center", fontsize=9.2, fontweight="bold", color="#0f172a", zorder=7)


def build_figure_2() -> Path:
    """Build 100% original, print-readable 3-row portrait schematic of the electro-metabolic oscillator."""
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(8.6, 11.0), dpi=240, facecolor="white")
    gs = gridspec.GridSpec(
        4, 1,
        height_ratios=[1.0, 1.0, 1.0, 0.50],
        hspace=0.10,
        left=0.02, right=0.98, top=0.935, bottom=0.015,
    )

    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])
    ax_bar = fig.add_subplot(gs[3])

    # ==================== PANEL A: BASAL POISED STATE ====================
    draw_wide_stage_frame(
        axA,
        bg_ext="#f0f9ff",
        stage_title="A   1. Basal Poised State (Brock 2001; Carr 2002)",
        badge_col="#0284c7",
        ext_caption="25–33 °C Rest  •  Isosmotic (~300 mOsm/kg)",
        readout_text="Vm ≈ −45 mV    |    [Ca²⁺]ᵢ ≈ 100 nM    |    pHᵢ ≈ 7.30    |    MICU1 Gate: CLOSED",
    )

    draw_channel_wide(axA, 12.5, 65.5, "#c4b5fd", "#6d28d9", open_pore=True, label="TRPM8")
    draw_channel_wide(axA, 29.0, 65.5, "#93c5fd", "#1d4ed8", open_pore=False, label="Naᵥ1.8", sublabel="Poised at\nactivation knee")
    draw_channel_wide(axA, 45.5, 65.5, "#fdba74", "#c2410c", open_pore=False, label="Caᵥ1.2 / 1.3", sublabel="Subthreshold /\nCaᵥ1.3 window")

    # Small tonic TRPM8 leak arrow
    axA.annotate("", xy=(12.5, 43.0), xytext=(12.5, 76.5),
                 arrowprops=dict(arrowstyle="->", color="#6d28d9", lw=2.0, ls="--"), zorder=7)
    axA.text(12.5, 40.5, "Small basal\nNa⁺/Ca²⁺ leak", ha="center", va="top", fontsize=8.6, fontweight="bold", color="#5b21b6")

    draw_mitochondrion_wide(axA, 76.5, 38.5, 17.5, 16.5, ca_loaded=False, gate_open=False, h_release=False)

    # ==================== PANEL B: OSMO-THERMAL BURST & L-TYPE AMPLIFICATION ====================
    draw_wide_stage_frame(
        axB,
        bg_ext="#e0f2fe",
        stage_title="B   2. Burst & Terminal L-Type Ca²⁺ (Gover 2003; Xu 2001)",
        badge_col="#7c3aed",
        ext_caption="Cooling (ΔT < 0) AND/OR +ΔOsm",
        readout_text="AP Bursts → −14 mV    |    [Ca²⁺]ᵢ: 100 → 515 nM    |    >300 nM MICU1 Set-Point Crossed",
    )

    draw_channel_wide(axB, 12.5, 65.5, "#a78bfa", "#5b21b6", open_pore=True, label="TRPM8")
    draw_channel_wide(axB, 29.0, 65.5, "#38bdf8", "#0369a1", open_pore=True, label="Naᵥ1.8")
    draw_channel_wide(axB, 45.5, 65.5, "#fb923c", "#c2410c", open_pore=True, label="Caᵥ1.2 / 1.3")

    # Strong ion influx arrows
    axB.annotate("", xy=(12.5, 40.0), xytext=(12.5, 76.5), arrowprops=dict(arrowstyle="->", color="#6d28d9", lw=2.8), zorder=7)
    axB.annotate("", xy=(29.0, 40.0), xytext=(29.0, 76.5), arrowprops=dict(arrowstyle="->", color="#0284c7", lw=2.8), zorder=7)
    axB.annotate("", xy=(45.5, 38.0), xytext=(45.5, 76.5), arrowprops=dict(arrowstyle="->", color="#ea580c", lw=3.2), zorder=7)

    # Horizontal depolarization coupling arrows
    axB.annotate("", xy=(24.5, 54.0), xytext=(16.5, 54.0), arrowprops=dict(arrowstyle="->", color="#0f172a", lw=1.9), zorder=8)
    axB.text(20.5, 49.5, "ΔVm", ha="center", va="top", fontsize=8.8, fontweight="bold", color="#0f172a")
    axB.annotate("", xy=(41.0, 54.0), xytext=(33.0, 54.0), arrowprops=dict(arrowstyle="->", color="#0f172a", lw=1.9), zorder=8)
    axB.text(37.0, 49.5, "AP spikes", ha="center", va="top", fontsize=8.4, fontweight="bold", color="#0f172a")

    # Ca2+ cloud rushing from Cav1.2/1.3 into open MICU1-MCU
    rng = np.random.default_rng(7)
    ca_x = rng.normal(53.5, 3.0, 14)
    ca_y = rng.normal(45.5, 4.8, 14)
    axB.scatter(ca_x, ca_y, s=30, color="#dc2626", edgecolor="white", lw=0.5, zorder=8)
    axB.annotate("", xy=(63.0, 55.0), xytext=(48.0, 45.0),
                 arrowprops=dict(arrowstyle="->", color="#dc2626", lw=2.6, connectionstyle="arc3,rad=-0.16"), zorder=9)
    axB.text(28.5, 29.0, "Terminal L-Type Ca²⁺ Surge (>300 nM; Scroggs 1992) → MICU1/MCU", ha="center", va="center",
             fontsize=8.2, fontweight="bold", color="#b91c1c",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#fef2f2", edgecolor="#fecaca", alpha=0.96), zorder=9)

    draw_mitochondrion_wide(axB, 76.5, 38.5, 17.5, 16.5, ca_loaded=True, gate_open=True, h_release=False)

    # ==================== PANEL C: MITOCHONDRIAL H+ RESET ====================
    draw_wide_stage_frame(
        axC,
        bg_ext="#fef2f2",
        stage_title="C   3. Mitochondrial H⁺ Reset (Werth 1994; Zhang 2025)",
        badge_col="#dc2626",
        ext_caption="H⁺ Shuts Cold/Osm TRPM8 (Spares Menthol!)",
        readout_text="Repolarizes → −51 mV    |    pHᵢ: 7.30 → 6.88 (TRPM8 Shut)    |    PMCA / NHE Clear → Next Burst",
    )

    draw_channel_wide(axC, 12.5, 65.5, "#ddd6fe", "#7c3aed", open_pore=False, label="TRPM8 (H⁺ Block)")
    draw_channel_wide(axC, 29.0, 65.5, "#bfdbfe", "#3b82f6", open_pore=False, label="Naᵥ1.8")
    draw_channel_wide(axC, 45.5, 65.5, "#fed7aa", "#c2410c", open_pore=False, label="Caᵥ1.2 / 1.3")

    draw_mitochondrion_wide(axC, 76.5, 38.5, 17.5, 16.5, ca_loaded=True, gate_open=True, h_release=True)

    # Cytosolic red H+ inhibition arc routing cleanly BELOW channels from TMBIM5 -> TRPM8
    axC.annotate("", xy=(14.5, 53.5), xytext=(78.0, 50.5),
                 arrowprops=dict(arrowstyle="-|>", color="#dc2626", lw=2.8, connectionstyle="arc3,rad=-0.20"), zorder=9)
    axC.plot([12.2, 16.8], [52.0, 55.2], color="#dc2626", lw=3.4, zorder=10)
    axC.text(35.0, 32.5, "Mitochondrial H⁺ Pulse (↑[H⁺]ᵢ, pHᵢ → 6.88)\nSelectively Closes Cold/Osmotic TRPM8 (Andersson 2004)",
             ha="center", va="center", fontsize=8.4, fontweight="bold", color="#b91c1c",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#fca5a5", alpha=0.97), zorder=10)

    # PMCA / NHE clearance arrow cleanly in right-hand extracellular margin below caption
    axC.annotate("", xy=(94.0, 72.5), xytext=(89.5, 57.5),
                 arrowprops=dict(arrowstyle="->", color="#0d9488", lw=2.3), zorder=8)
    axC.text(97.0, 75.0, "PMCA (Gover 2007)\n& NHE clear Ca²⁺/H⁺", ha="right", va="bottom",
             fontsize=7.6, fontweight="bold", color="#0f766e",
             bbox=dict(boxstyle="round,pad=0.18", facecolor="#f0fdfa", edgecolor="#99f6e4", alpha=0.96), zorder=9)

    # ==================== BOTTOM STRIP: 2x2 CLINICAL / GENETIC FAILURE NODES ====================
    ax_bar.set_xlim(0, 100)
    ax_bar.set_ylim(0, 100)
    ax_bar.axis("off")
    ax_bar.add_patch(FancyBboxPatch((0.6, 2.0), 98.8, 96.0, boxstyle="round,pad=0.4",
                                    facecolor="#f8fafc", edgecolor="#94a3b8", lw=1.2))
    ax_bar.text(50.0, 91.5, "Human Genetic, Metabolic & Pharmacological Failure Nodes (Figures 4–5)",
                ha="center", va="center", fontsize=9.2, fontweight="bold", color="#0f172a")

    pills = [
        (25.5, 65.0, "#ede9fe", "#5b21b6", "TRPM8 GoF (D665N / R30Q)\nEctopic LASIK & Trigeminal Neuralgia Bursting"),
        (74.5, 65.0, "#dbeafe", "#1d4ed8", "SCN10A GoF (I206M / G1662S)\nTonic Plateau Lock-In & Small-Fiber Neuropathy"),
        (25.5, 24.5, "#ffedd5", "#c2410c", "Caᵥ1.2 / 1.3 Dihydropyridine Pharmacology\nBay K 8644 / Nifedipine (+Vascular Confound)"),
        (74.5, 24.5, "#fee2e2", "#b91c1c", "MICU1 LoF (MPXPS) & Diabetic MCU KO\nMatrix Ca²⁺ Overload vs. Axonal Rescue"),
    ]
    for px, py, fc, ec, txt in pills:
        ax_bar.add_patch(FancyBboxPatch((px - 23.4, py - 15.5), 46.8, 31.0, boxstyle="round,pad=0.35",
                                        facecolor=fc, edgecolor=ec, lw=1.1))
        ax_bar.text(px, py, txt, ha="center", va="center", fontsize=8.2, fontweight="bold", color=ec, linespacing=1.25)

    fig.suptitle(
        "Figure 2. Electro-Metabolic Relaxation Oscillator in Cold/Osmo-Sensitive Terminals",
        fontsize=11.4,
        fontweight="bold",
        color="#0f172a",
        y=0.974,
    )

    out_path = OUT_DIR / "fig2_mechanistic_circuit.png"
    fig.savefig(out_path, dpi=240, facecolor="white")
    plt.close(fig)
    return out_path


def _render_eq_block(lines: list[tuple[str, str]], out_path: Path, figsize: tuple[float, float]) -> Path:
    """Render clean 300-DPI STIX display equations with right-aligned equation numbers."""
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["font.family"] = "STIXGeneral"
    fig, ax = plt.subplots(figsize=figsize, dpi=300, facecolor="white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    n = len(lines)
    ys = np.linspace(0.78, 0.22, n) if n > 1 else [0.5]
    for y, (eq, num) in zip(ys, lines):
        ax.text(0.47, y, eq, ha="center", va="center", fontsize=14.5, color="#0f172a")
        ax.text(0.96, y, num, ha="right", va="center", fontsize=12.5, color="#475569")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.02)
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    return out_path


def build_equation_blocks() -> list[Path]:
    """Build 300-DPI STIX display-equation blocks (Eqs. 1–7) for Google Docs import."""
    p_eq1 = _render_eq_block(
        [
            (r"$P_{\mathrm{open}}(V,\ T,\ \mathrm{Osm},\ \mathrm{pH}_i)\ =\ \left[\,1 + \exp\!\left(\dfrac{\Delta G_{\mathrm{total}}}{RT}\right)\right]^{-1}$", "(1)"),
            (r"$\Delta G_{\mathrm{total}}\ =\ \Delta H(T) - T\Delta S(T) - zFV + \Delta G_{\mathrm{osm}}(\mathrm{Osm}_{\mathrm{ext}} - \mathrm{Osm}_{\mathrm{int}}) + \Delta G_{\mathrm{tens}}(\gamma) + \Delta G_{\mathrm{H}^+}(\mathrm{pH}_i)$", "(2)"),
        ],
        OUT_DIR / "eq_block1_thermodynamics.png",
        figsize=(10.2, 1.55),
    )
    p_eq2 = _render_eq_block(
        [
            (r"$c_{\mathrm{osm}}\ =\ -\dfrac{\ln(a_w)}{M_w}, \qquad \Pi\ =\ -\dfrac{RT}{\bar{V}_w}\ln(a_w) \qquad \left(M_w = 0.018015\ \mathrm{kg\,mol}^{-1}\right)$", "(3)"),
            (r"$\Delta c_{\mathrm{osm,SG}}\ \approx\ -\dfrac{\Delta a_{w,\mathrm{SG}}}{M_w\,a_{w,0}}\ =\ -\dfrac{(1 - S_{\mathrm{SC}})\,\Delta a_{w,\mathrm{surf}}}{M_w\,a_{w,0}}, \qquad \left(\dfrac{dc_{\mathrm{osm}}}{da_w}\right)_{\!a_{w,0}} =\ -55.8\ \mathrm{Osm\,kg}^{-1}$", "(4)"),
        ],
        OUT_DIR / "eq_block2_water_activity.png",
        figsize=(10.2, 1.65),
    )
    p_eq3 = _render_eq_block(
        [
            (r"$C_m\,\dfrac{dV}{dt}\ =\ -k\,\left[\,I_{\mathrm{TRPM8}}(V, T, \mathrm{Osm}, s_{\mathrm{pH}}) + I_{\mathrm{Na_v1.8}}(V) + I_{\mathrm{Ca_L}}(V, c) + I_L(V, c)\,\right]$", "(5)"),
            (r"$\dfrac{dc}{dt}\ =\ k\,\left[\,J_{\mathrm{in,M8}} + J_{\mathrm{in,CaL}} - J_{\mathrm{PMCA}}(c) - J_{\mathrm{MCU}}(c, \Psi_m) + J_{\mathrm{m,out}}(c_m)\,\right]$", "(6)"),
            (r"$\dfrac{d(\Delta[\mathrm{H}^+]_i)}{dt}\ =\ k\left(\alpha_H J_{\mathrm{MCU}}(c, \Psi_m) - \dfrac{\Delta[\mathrm{H}^+]_i}{\tau_H}\right), \qquad \dfrac{ds_{\mathrm{pH}}}{dt}\ =\ k\,\dfrac{s_\infty([\mathrm{H}^+]_i) - s_{\mathrm{pH}}}{\tau_s}$", "(7)"),
        ],
        OUT_DIR / "eq_block3_ode_system.png",
        figsize=(10.2, 2.35),
    )
    # Also keep eq_block2_ode_system.png synced as an alias in case any legacy link references it
    shutil.copy2(p_eq3, OUT_DIR / "eq_block2_ode_system.png")
    return [p_eq1, p_eq2, p_eq3]


def main():
    p1 = build_figure_1()
    p2 = build_figure_2()
    eqs = build_equation_blocks()
    print("Built copyright-clean Figure 1:", p1)
    print("Built copyright-clean Figure 2:", p2)
    for ep in eqs:
        print("Built STIX equation block:", ep)


if __name__ == "__main__":
    main()
