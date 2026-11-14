#!/usr/bin/env python3
"""Generate two new figures:
1. fig_variants.pdf — distinction between t_4^{(5)}, t_4, T_4 (for §1.1)
2. fig_T4_60.pdf — actual geometric N=60 solution (replacing param-space heatmap)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import json

# ============================================================
# FIGURE 1: Three definitions distinction
# ============================================================
def make_variants():
    fig, axes = plt.subplots(1, 3, figsize=(9, 2.8), sharex=False, sharey=False)
    fig.subplots_adjust(wLeft := 0.05, right=0.98, wspace=0.18, bottom=0.05, top=0.88)

    # Common geometry: 9 points arranged to give one 5-line (horizontal y=0) and one 4-line (diagonal)
    # Points for 5-line: 5 points on y=0, x=0..4
    # Points for 4-line: 4 points on y=x+offset, but one overlaps with the 5-line? Let's design so
    # we have clean distinction: 9 points total, with:
    #   H: (0,0),(1,0),(2,0),(3,0),(4,0)  -> 5 collinear (horizontal)
    #   D: (0,1),(1,2),(2,3),(4,0) last point shared -> Wait need distinct counts
    # Simpler: total 10 points, with a 5-line and a 4-line that share no points?
    # Let's do: 9 points = 5 on horizontal + 4 on diagonal, but arrange so diagonal has 4 points
    # including one that is NOT on horizontal, but we can illustrate sharing with 0 overlap for clarity.

    # Define points
    H = [(i, 0) for i in range(5)]  # horizontal line y=0, 5 points
    D = [(i+0.3, 1.2+i*0.7) for i in range(4)]  # diagonal line, 4 points, offset above
    # plus one isolated point to make it 10? Actually H+D = 9 points, fine
    all_pts = H + D

    # For each panel, we draw same geometry but highlight differently
    configs = [
        {
            "title": r"$t_4^{(5)}(n)$: no $5$ collinear",
            "subtitle": "5-line forbidden, not counted",
            "highlight": "forbidden",  # show 5-line as forbidden
            "count": 1,
            "count_label": "counts 1 four-line",
        },
        {
            "title": r"$t_4(n)$: count exactly $4$",
            "subtitle": "5-line exists but not counted",
            "highlight": "ignore",
            "count": 1,
            "count_label": "counts 1 four-line",
        },
        {
            "title": r"$T_4(n)$: count $\geq 4$",
            "subtitle": "both lines counted",
            "highlight": "both",
            "count": 2,
            "count_label": "counts 2 lines",
        },
    ]

    for ax, cfg in zip(axes, configs):
        # Draw faint grid
        ax.set_xlim(-0.7, 5.0)
        ax.set_ylim(-0.8, 4.5)
        ax.set_aspect('equal')
        ax.axis('off')

        # Draw the 5-point horizontal line
        if cfg["highlight"] == "forbidden":
            # red dashed, with X
            ax.plot([-0.5, 4.7], [0, 0], color='#e74c3c', linewidth=2.2, linestyle='--', zorder=1, alpha=0.9)
            # red X at middle
            ax.plot([2-0.25, 2+0.25], [0-0.25, 0+0.25], color='#e74c3c', linewidth=2, zorder=4)
            ax.plot([2-0.25, 2+0.25], [0+0.25, 0-0.25], color='#e74c3c', linewidth=2, zorder=4)
            ax.text(2, -0.55, "5 points — forbidden", ha='center', va='top', fontsize=6.5, color='#c0392b', style='italic')
        elif cfg["highlight"] == "ignore":
            # gray, ignored
            ax.plot([-0.5, 4.7], [0, 0], color='#bdc3c7', linewidth=1.8, linestyle=':', zorder=1)
            ax.text(2, -0.55, "5 points — not counted", ha='center', va='top', fontsize=6.5, color='#7f8c8d', style='italic')
        else:  # both
            ax.plot([-0.5, 4.7], [0, 0], color='#2980b9', linewidth=2.2, zorder=1)
            ax.text(2, -0.55, "5 points — counted", ha='center', va='top', fontsize=6.5, color='#2471a3', weight='bold')

        # Draw the 4-point diagonal line
        # diagonal line coords
        d_x = [p[0] for p in D]
        d_y = [p[1] for p in D]
        # extend slightly
        ax.plot([d_x[0]-0.3, d_x[-1]+0.3], [d_y[0]-0.2, d_y[-1]+0.2], color='#27ae60', linewidth=2.2, zorder=1)
        ax.text(np.mean(d_x), np.mean(d_y)+0.55, "4 points — counted", ha='center', va='bottom', fontsize=6.5, color='#1e8449', weight='bold')

        # Draw points
        for (x, y) in H:
            ax.plot(x, y, 'o', markersize=7, markerfacecolor='white', markeredgecolor='#2c3e50', markeredgewidth=1.2, zorder=3)
        for (x, y) in D:
            ax.plot(x, y, 'o', markersize=7, markerfacecolor='white', markeredgecolor='#2c3e50', markeredgewidth=1.2, zorder=3)

        # Count badge at top
        ax.text(0.5, 0.96, f"= {cfg['count']}", transform=ax.transAxes, ha='center', va='top',
                fontsize=13, weight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='#fef9e7' if cfg["count"]==1 else '#eaf2f8', edgecolor='#aab7b8', linewidth=0.8))

        ax.set_title(cfg["title"] + "\n" + r"$\,$" + cfg["subtitle"], fontsize=7.5, pad=8, color='#2c3e50')

    # Overall title
    fig.suptitle(r"Three ways to count: the role of a $5$-point line", fontsize=9, weight='bold', y=0.99, color='#1a1a1a')
    # Add legend at bottom
    fig.text(0.5, 0.01, r"All three coincide when no line has $\geq 5$ points; they differ precisely when long rows are present.  The constructions of \S 5 use $t_4^{(5)}$, while \S 6 uses $T_4$.",
             ha='center', fontsize=6.2, color='#666666', style='italic')

    fig.savefig('paper/figures/fig_variants.pdf', bbox_inches='tight', pad_inches=0.08)
    fig.savefig('paper/figures/fig_variants.png', bbox_inches='tight', pad_inches=0.08, dpi=200)
    print("wrote fig_variants.pdf/png")

# ============================================================
# FIGURE 2: Actual N=60 geometric solution
# ============================================================
def make_T4_60_geometric():
    # Load certificate for A,B,C,D intervals
    d = json.load(open('atleast4_N60_certificate.json'))
    A, B, C, D = d['A'], d['B'], d['C'], d['D']
    # Generate theoretical points in natural coordinates (before scaling)
    # Points: (1/(180+x), 0), (0, 1/(60+y)), (1/(240+z), 1/(240+z)), (1/(120+w), -1/(120+w))
    pts_A = [(1/(180+x), 0) for x in A]
    pts_B = [(0, 1/(60+y)) for y in B]
    pts_C = [(1/(240+z), 1/(240+z)) for z in C]
    pts_D = [(1/(120+w), -1/(120+w)) for w in D]

    # For plotting, scale for visibility: multiply by 30000 to spread
    # Actually 1/180 ~0.0055, 1/60~0.016, so range ~0.003-0.016, scale 3000 gives 10-50 units
    scale = 6000  # tuned for nice spread

    def sc(p): return (p[0]*scale, p[1]*scale)

    pts_A_s = [sc(p) for p in pts_A]
    pts_B_s = [sc(p) for p in pts_B]
    pts_C_s = [sc(p) for p in pts_C]
    pts_D_s = [sc(p) for p in pts_D]

    all_pts = pts_A_s + pts_B_s + pts_C_s + pts_D_s
    n = len(all_pts)
    assert n == 60, n

    # Build lookup for transversals: condition x+y in C and x-y in D where x in A, y in B
    # A,B here are denominators 180+x, 60+y etc - but we have original interval values A,B
    # Let's use original integer parameters: A_vals = 0..12 mapping to 180+x, etc
    # Condition is (x_param + y_param) in C_params and (x_param - y_param) in D_params
    # where x_param = x offset (0..12), y_param = y offset (0..12)
    # C = 4..20, D = -8..8
    A_off = list(range(13))  # 0..12
    B_off = list(range(13))
    C_set = set(range(4,21))  # 4..20
    D_set = set(range(-8,9))  # -8..8
    # Build all 4-point transversals as point indices
    # Each transversal uses 1 point from A, 1 from B, 1 from C, 1 from D? No - need to verify
    # The 4-axis construction: a transversal with parameters (x,y) meets axis C at p3x+q3y and D at p4x+q4y
    # For original slopes (1,1) and (1,-1), C param = x+y, D param = x-y
    # So for each (x,y) that satisfies, we get 4 points: A[x], B[y], C[x+y], D[x-y]
    # That's a 4-point line.

    # Collect lines as tuples of 4 point indices (in the global all_pts list)
    # A indices 0..12, B 13..25, C 26..42, D 43..59
    lines = []
    for xi, x in enumerate(A_off):
        for yi, y in enumerate(B_off):
            if (x+y) in C_set and (x-y) in D_set:
                # map to C index: C holds values 244..260 = 240+(z) where z in 4..20
                # C values are denominators, but our C_set is 4..20 (the z values)
                # C index = (x+y) - 4
                ci = (x+y) - 4  # 0..16
                # D index: D values -8..8 map to D list 112..128 = 120+w
                di = (x-y) + 8  # 0..16
                lines.append((xi, 13+yi, 26+ci, 43+di))
    print(f"found {len(lines)} 4-point transversals (expected 129)")

    # Now plot
    fig, ax = plt.subplots(figsize=(6.0, 5.2), dpi=200)
    ax.set_aspect('equal')

    # Draw axes faintly - four lines through origin
    # x-axis (y=0), y-axis (x=0), y=x, y=-x
    # Extend beyond points
    lim = 38
    ax.plot([-2, lim], [0, 0], color='#bdc3c7', linewidth=0.8, linestyle='-', zorder=0)
    ax.plot([0, 0], [-lim, lim], color='#bdc3c7', linewidth=0.8, linestyle='-', zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color='#bdc3c7', linewidth=0.8, linestyle='-', zorder=0)
    ax.plot([-lim, lim], [lim, -lim], color='#bdc3c7', linewidth=0.8, linestyle='-', zorder=0)

    # Draw all 4-point transversal lines faintly (light blue, thin)
    for (ia, ib, ic, id_) in lines:
        pts_line = [all_pts[ia], all_pts[ib], all_pts[ic], all_pts[id_]]
        # Sort by x for line drawing (collinear, but just draw segment between extremes)
        xs = [p[0] for p in pts_line]
        ys = [p[1] for p in pts_line]
        # Draw as thin line
        ax.plot(xs, ys, color='#aed6f1', linewidth=0.35, alpha=0.45, zorder=1, solid_capstyle='round')

    # Draw points colored by axis
    colors = {
        'A': '#e74c3c',  # red - x-axis
        'B': '#2980b9',  # blue - y-axis
        'C': '#27ae60',  # green - diagonal
        'D': '#8e44ad',  # purple - anti-diagonal
    }
    # A: red squares
    for p in pts_A_s:
        ax.plot(p[0], p[1], 's', markersize=4.5, markerfacecolor=colors['A'], markeredgecolor='white', markeredgewidth=0.6, zorder=5)
    for p in pts_B_s:
        ax.plot(p[0], p[1], 'o', markersize=4.5, markerfacecolor=colors['B'], markeredgecolor='white', markeredgewidth=0.6, zorder=5)
    for p in pts_C_s:
        ax.plot(p[0], p[1], 'D', markersize=4.5, markerfacecolor=colors['C'], markeredgecolor='white', markeredgewidth=0.6, zorder=5)
    for p in pts_D_s:
        ax.plot(p[0], p[1], '^', markersize=5, markerfacecolor=colors['D'], markeredgecolor='white', markeredgewidth=0.6, zorder=5)

    # Highlight a couple of example 4-point lines in bold (pick 2 that are well separated)
    example_indices = [0, 40, 80]  # pick a few
    highlight_colors = ['#c0392b', '#1a5276', '#196f3d']
    for idx, hi in zip(example_indices, highlight_colors):
        if idx < len(lines):
            ia, ib, ic, id_ = lines[idx]
            pts_line = [all_pts[ia], all_pts[ib], all_pts[ic], all_pts[id_]]
            xs = [p[0] for p in pts_line]
            ys = [p[1] for p in pts_line]
            ax.plot(xs, ys, color=hi, linewidth=1.4, alpha=0.9, zorder=4)
            # mark those 4 points with larger highlight
            for (x, y) in pts_line:
                ax.plot(x, y, 'o', markersize=6, markerfacecolor='none', markeredgecolor=hi, markeredgewidth=1.3, zorder=6)

    ax.set_xlim(-3, 38)
    ax.set_ylim(-34, 38)
    ax.set_xlabel(r'$x$ (scaled)', fontsize=8)
    ax.set_ylabel(r'$y$ (scaled)', fontsize=8)
    ax.tick_params(labelsize=7, length=2.5)

    # Title and subtitle
    ax.set_title(r"$T_4(60)\geq133$: 60 points on four concurrent axes and 129 transversals", fontsize=9, weight='bold', pad=10)
    # Legend
    legend_elements = [
        Line2D([0],[0], marker='s', color='w', markerfacecolor=colors['A'], markersize=6, label=r'$A$: $y=0$ (13 pts)', markeredgecolor='#888', markeredgewidth=0.5),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=colors['B'], markersize=6, label=r'$B$: $x=0$ (13 pts)', markeredgecolor='#888', markeredgewidth=0.5),
        Line2D([0],[0], marker='D', color='w', markerfacecolor=colors['C'], markersize=6, label=r'$C$: $y=x$ (17 pts)', markeredgecolor='#888', markeredgewidth=0.5),
        Line2D([0],[0], marker='^', color='w', markerfacecolor=colors['D'], markersize=6, label=r'$D$: $y=-x$ (17 pts)', markeredgecolor='#888', markeredgewidth=0.5),
    ]
    leg = ax.legend(handles=legend_elements, loc='upper left', fontsize=6.5, frameon=True, facecolor='white', edgecolor='#bdc3c7', borderpad=0.4, handletextpad=0.4)
    leg.get_frame().set_linewidth(0.6)

    # Annotation box
    props = dict(boxstyle="round,pad=0.35", facecolor='#fef9e7', edgecolor='#d4ac0d', linewidth=0.7)
    ax.text(0.98, 0.02, "Each pale blue line is a 4-point transversal\n"
                        "($x+y\\in C$, $x-y\\in D$): 129 in total.\n"
                        "Bold lines highlight three examples.\n"
                        "Two axes carry 17 points ($C,D$) → not admissible\n"
                        "for $t_4^{(5)}$ (counts toward $T_4$).",
             transform=ax.transAxes, ha='right', va='bottom', fontsize=6.2,
             bbox=props, color='#2c3e50', linespacing=1.3)

    fig.tight_layout(pad=0.6)
    fig.savefig('paper/figures/fig_T4_60.pdf', bbox_inches='tight', pad_inches=0.15)
    fig.savefig('paper/figures/fig_T4_60.png', bbox_inches='tight', pad_inches=0.15, dpi=200)
    print("wrote fig_T4_60.pdf/png (geometric)")

if __name__ == '__main__':
    make_variants()
    make_T4_60_geometric()
