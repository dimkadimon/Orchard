#!/usr/bin/env python3
"""Fix fig_variants only: new notation s_4 and visible top text."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def make_variants():
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.0), sharex=False, sharey=False)
    # more room on top and bottom
    fig.subplots_adjust(left=0.05, right=0.98, wspace=0.22, bottom=0.18, top=0.82)

    H = [(i, 0) for i in range(5)]
    D = [(i+0.3, 1.2+i*0.7) for i in range(4)]

    configs = [
        {
            "title": r"$s_4(n)$: strict, no $5$ collinear",
            "subtitle": "5-line forbidden, not counted",
            "highlight": "forbidden",
            "count": 1,
        },
        {
            "title": r"$t_4(n)$: count exactly $4$",
            "subtitle": "5-line exists but not counted",
            "highlight": "ignore",
            "count": 1,
        },
        {
            "title": r"$T_4(n)$: count $\geq 4$",
            "subtitle": "both lines counted",
            "highlight": "both",
            "count": 2,
        },
    ]

    for ax, cfg in zip(axes, configs):
        ax.set_xlim(-0.7, 5.0)
        ax.set_ylim(-0.8, 4.5)
        ax.set_aspect('equal')
        ax.axis('off')

        if cfg["highlight"] == "forbidden":
            ax.plot([-0.5, 4.7], [0, 0], color='#e74c3c', linewidth=2.2, linestyle='--', zorder=1, alpha=0.9)
            ax.plot([2-0.25, 2+0.25], [0-0.25, 0+0.25], color='#e74c3c', linewidth=2, zorder=4)
            ax.plot([2-0.25, 2+0.25], [0+0.25, 0-0.25], color='#e74c3c', linewidth=2, zorder=4)
            ax.text(2, -0.55, "5 points — forbidden", ha='center', va='top', fontsize=6.5, color='#c0392b', style='italic')
        elif cfg["highlight"] == "ignore":
            ax.plot([-0.5, 4.7], [0, 0], color='#bdc3c7', linewidth=1.8, linestyle=':', zorder=1)
            ax.text(2, -0.55, "5 points — not counted", ha='center', va='top', fontsize=6.5, color='#7f8c8d', style='italic')
        else:
            ax.plot([-0.5, 4.7], [0, 0], color='#2980b9', linewidth=2.2, zorder=1)
            ax.text(2, -0.55, "5 points — counted", ha='center', va='top', fontsize=6.5, color='#2471a3', weight='bold')

        d_x = [p[0] for p in D]
        d_y = [p[1] for p in D]
        ax.plot([d_x[0]-0.3, d_x[-1]+0.3], [d_y[0]-0.2, d_y[-1]+0.2], color='#27ae60', linewidth=2.2, zorder=1)
        ax.text(np.mean(d_x), np.mean(d_y)+0.55, "4 points — counted", ha='center', va='bottom', fontsize=6.5, color='#1e8449', weight='bold')

        for (x, y) in H:
            ax.plot(x, y, 'o', markersize=7, markerfacecolor='white', markeredgecolor='#2c3e50', markeredgewidth=1.2, zorder=3)
        for (x, y) in D:
            ax.plot(x, y, 'o', markersize=7, markerfacecolor='white', markeredgecolor='#2c3e50', markeredgewidth=1.2, zorder=3)

        ax.text(0.5, 0.96, f"= {cfg['count']}", transform=ax.transAxes, ha='center', va='top',
                fontsize=13, weight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='#fef9e7' if cfg["count"]==1 else '#eaf2f8', edgecolor='#aab7b8', linewidth=0.8))

        ax.set_title(cfg["title"] + "\n" + r"$\,$" + cfg["subtitle"], fontsize=7.8, pad=6, color='#2c3e50')

    fig.suptitle(r"Three ways to count: the role of a $5$-point line", fontsize=9.5, weight='bold', y=0.96, color='#1a1a1a')
    fig.text(0.5, 0.05, r"All three coincide when no line has $\geq 5$ points; they differ precisely when long rows are present.  The constructions of \S 5 use $s_4$, while \S 6 uses $T_4$.",
             ha='center', fontsize=6.4, color='#666666', style='italic', wrap=True)

    fig.savefig('paper/figures/fig_variants.pdf', bbox_inches='tight', pad_inches=0.08)
    fig.savefig('paper/figures/fig_variants.png', bbox_inches='tight', pad_inches=0.08, dpi=200)
    print("wrote fig_variants.pdf/png with fixed top and s_4 notation")

if __name__ == '__main__':
    make_variants()
