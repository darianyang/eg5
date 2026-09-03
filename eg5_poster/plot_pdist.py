
import wedap
import matplotlib.pyplot as plt

plt.style.use("~/Apps/wedap/wedap/styles/default.mplstyle")

def plotter(system="wt", h5="west.h5"):

    plot_options = {"h5": f"eg5_{system}/{h5}", 
                    "data_type": "average", 
                    "plot_mode": "hexbin3d",
                    "Xname": "pcoord",
                    "Xindex": 0,
                    "xlim": (-5, 60),
                    "xlabel": r"ADP RMSD ($\AA$)",
                    "Yname": "pcoord",
                    "Yindex": 1,
                    "ylim": (-250, 410),
                    "ylabel": r"Interaction Energy ($\frac{kcal}{mol}$)",
                    "Zname": "pcoord",
                    "Zindex": 2,
                    "clim": (2, 10),
                    "cbar_label": r"Min Distance ($\AA$)",
                    "cmap": "bwr",
                    }

    wp = wedap.H5_Plot(**plot_options).plot()
    plt.savefig(f"pdist_{system}_{h5.replace('.h5', '')}.png", dpi=300, bbox_inches="tight")

plotter("wt", h5="west.h5")
plotter("wt", h5="succ_only.h5")
plotter("mon", h5="west.h5")
plotter("mon", h5="succ_only.h5")