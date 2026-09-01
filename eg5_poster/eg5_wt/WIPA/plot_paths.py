import wedap
import matplotlib.pyplot as plt

plt.style.use("~/github/wedap/wedap/styles/default.mplstyle")

clusters = [(548, 302), (461, 293), (468, 6), (526, 300)]

wedap_plot = wedap.H5_Plot(h5="succ_only.h5", data_type="average", Xname="pcoord", Xindex=0, Yname="pcoord", Yindex=1, Zname="pcoord", Zindex=2, plot_mode="hexbin3d", p_max=12, p_min=2, cmap="copper")#, cbar_label="ADP RMSD ($\AA$)")

wedap_plot.plot()
#wedap_plot.ax.set(xlim=(-1,55), ylim=(-250, 420), xlabel=r"Interaction Energy ($\frac{kcal}{mol}$)", ylabel="Min Distance ($\AA$)")
colors = ["#377eb8", "#d62728", "#4daf4a", "#f781bf"]
for i, c in enumerate(clusters):
    wedap_plot.plot_trace(c, ax=wedap_plot.ax, linewidth=0.5, color=colors[i])

plt.tight_layout()
plt.savefig("pathways_traced2.pdf")
plt.show()