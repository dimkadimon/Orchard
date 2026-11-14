"""Diagrams for k>4 improvements: classical vs second-order upper bound."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from math import comb
from second_order_k import melchior_k, bound

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':0.7})

colors_classical = '#bdc3c7'
colors_new = {'5':'#2980b9','6':'#27ae60','7':'#8e44ad','8':'#e67e22'}
ks = [5,6,7,8]

# Create 2x2 grid
fig, axes = plt.subplots(2,2, figsize=(8.2, 6.0), sharex=False, sharey=False)
fig.subplots_adjust(left=0.07, right=0.98, wspace=0.22, hspace=0.45, top=0.90, bottom=0.10)
axes = axes.flatten()

for idx,k in enumerate(ks):
    ax = axes[idx]
    ns = list(range(k+1, 80 if k<=6 else 85))
    classical = [melchior_k(n,k) for n in ns]
    new = [bound(n,k) for n in ns]
    # Plot
    ax.plot(ns, classical, color=colors_classical, linewidth=1.6, linestyle='--', label='classical (Prop. 4.1)', zorder=1)
    ax.plot(ns, new, color=colors_new[str(k)], linewidth=2.0, label='new (Cor. 4.3)', zorder=2)
    # Shade improvement
    ax.fill_between(ns, new, classical, where=[c>n for c,n in zip(classical,new)], color=colors_new[str(k)], alpha=0.18, zorder=0)
    ax.set_title(f'$k={k}$', fontsize=10, weight='bold', pad=6)
    ax.set_xlabel('$n$', fontsize=8)
    ax.set_ylabel(f'$T_{k}(n)$ upper bound', fontsize=8)
    ax.set_xlim(min(ns), max(ns))
    ax.set_ylim(0, max(classical)*1.08)
    ax.grid(True, color='#ecf0f1', linewidth=0.5, zorder=-1)
    ax.tick_params(labelsize=7)
    # Annotate max gain
    max_gap = max(c-n for c,n in zip(classical,new))
    n_max = ns[[c-n for c,n in zip(classical,new)].index(max_gap)]
    if max_gap>0:
        ax.annotate(f'max gain {max_gap}', xy=(n_max, new[ns.index(n_max)]), xytext=(8,12), textcoords='offset points',
                    ha='center', fontsize=6.5, color=colors_new[str(k)], weight='bold',
                    bbox=dict(boxstyle="round,pad=0.2", fc='white', ec=colors_new[str(k)], alpha=0.9, linewidth=0.6),
                    arrowprops=dict(arrowstyle='->', color=colors_new[str(k)], lw=0.8))
    leg = ax.legend(fontsize=6.5, loc='upper left', frameon=True, facecolor='white', edgecolor='#bdc3c7', handlelength=1.8)
    leg.get_frame().set_linewidth(0.6)

fig.suptitle(r'Second-order bound strictly improves the classical Melchior bound for $k\geq5$' '\n' r'Classical $\, \lfloor(\binom{n}{2}-3)/(\binom{k}{2}+k-3)\rfloor$ vs new (Cor.~\ref{cor:test}), $n$ up to 80',
             fontsize=10, weight='bold', y=0.97)

# add table inset summary
fig.text(0.5, 0.01,
         r'Table~\ref{tab:k}: improves to $n=25,37,51,67$ for $k=5,6,7,8$. Largest early gain: $k=8$ at $n=50$: $37\to19$ ($-49\%$).',
         ha='center', fontsize=7, color='#555555', style='italic')

fig.savefig('paper/figures/fig_k_bounds.pdf', bbox_inches='tight', pad_inches=0.1)
fig.savefig('paper/figures/fig_k_bounds.png', bbox_inches='tight', pad_inches=0.1, dpi=200)
print("wrote fig_k_bounds.pdf/png")

# Single combined overlay for talk
fig2, ax2 = plt.subplots(figsize=(6.5, 4.0))
fig2.subplots_adjust(left=0.10, right=0.96, bottom=0.12, top=0.88)
for k in ks:
    ns = list(range(k+1, 75))
    classical = [melchior_k(n,k) for n in ns]
    new = [bound(n,k) for n in ns]
    ax2.plot(ns, [c-n for c,n in zip(classical,new)], color=colors_new[str(k)], linewidth=1.9, label=f'$k={k}$')
ax2.set_xlabel('$n$', fontsize=9)
ax2.set_ylabel('improvement $\\Delta = $ classical $-$ new', fontsize=9)
ax2.set_title('Gain of second-order test over classical bound', fontsize=10, weight='bold')
ax2.grid(True, color='#ecf0f1', linewidth=0.5)
ax2.legend(fontsize=8)
fig2.savefig('paper/figures/fig_k_gain.pdf', bbox_inches='tight', pad_inches=0.1)
fig2.savefig('paper/figures/fig_k_gain.png', bbox_inches='tight', pad_inches=0.1, dpi=200)
print("wrote fig_k_gain.pdf/png")
