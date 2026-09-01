"""
Plot the fluxes and rates from direct.h5 files.

list(h5):
['arrivals', 'avg_color_probs', 'avg_conditional_fluxes', 'avg_rates', 
'avg_state_probs', 'avg_total_fluxes', 'color_prob_evolution', 'conditional_arrivals',
'conditional_flux_evolution', 'conditional_fluxes', 'duration_count', 'durations', 
'rate_evolution', 'state_labels', 'state_pop_evolution', 'target_flux_evolution', 'total_fluxes']

  /target_flux_evolution [window,state]
    Total flux into a given macro state based on
    windows of iterations of varying width, as in /rate_evolution.
  /conditional_flux_evolution [window,state,state]
    State-to-state fluxes based on windows of
    varying width, as in /rate_evolution.

The structure of these datasets is as follows:
  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).
  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).
  expected
    (Floating-point) Expected (mean) value of the observable as evaluated within
    this window, in units of inverse tau.
  ci_lbound
    (Floating-point) Lower bound of the confidence interval of the observable
    within this window, in units of inverse tau.
  ci_ubound
    (Floating-point) Upper bound of the confidence interval of the observable
    within this window, in units of inverse tau.
  stderr
    (Floating-point) The standard error of the mean of the observable
    within this window, in units of inverse tau.
  corr_len
    (Integer) Correlation length of the observable within this window, in units
    of tau.
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py

#from scipy.signal import savgol_filter

plt.style.use("~/github/wedap/wedap/styles/default.mplstyle")

# TODO: function to make 4 panel plot
    # Plot of P_A, P_B, rate_AB, rate_BA, all as function of WE iteration

class Kinetics:

    def __init__(self, min_state=None, max_state=None, scheme=None, prefix=None,
                 tau=10**-10, state=1, label=None, units="rates", ax=None, savefig=False):
        """ TODO
        Parameters
        ----------
        scheme : str
            Scheme/directory for the direct.h5 file.
        tau : float
            The resampling interval of the WE simualtion.
            This should be in seconds, default 100ps = 100 * 10^-12 = 10^-10
        state : int
            State for flux calculations, 0 = A and 1 = B.
        label : str
            Data label.
        units : str
            Can be `rates` (default) or `mfpts`.
        ax : mpl axes object
        savefig : bool
            Wether or not to save the figure to a png.
        """
        self.scheme = scheme
        self.prefix = prefix
        self.tau = tau
        self.state = state
        self.label = label
        self.units = units

        self.min_state = min_state
        self.max_state = max_state
        if min_state and max_state:
            self.states = [i for i in range(min_state, max_state, 1)]
            self.n_states = len(self.states)
        else:
            self.states = None
            self.n_states = 1

        # create new fig or plot onto existing
        if ax is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = ax
            self.fig = plt.gcf()

        self.savefig = savefig

    def extract_rate(self):
        """
        Get the raw rate array from one direct.h5 file.
        """
        print(self.scheme)
        # read in direct.h5 file
        h5 = h5py.File(f"{self.scheme}/direct.h5", "r")

        # flux evolution dataset from cumulative evolution mode:
        # When calculating time evolution of rate estimates, 
        # ``cumulative`` evaluates rates over windows starting with --start-iter and 
        # getting progressively wider to --stop-iter by steps of --step-iter.
        fluxes = np.array(h5["target_flux_evolution"])

        # conditional fluxes are macrostate to macrostate
        # 2 dimensions: [(0 -> 0, 0 -> 1), 
        #                (1 -> 0, 1 -> 1)] 
        # I want 0 -> 1
        #fluxes = np.array(h5["conditional_flux_evolution"])[:,:,1]

        # third column (expected) of the state (A/0 or B/1) flux dataset (flux into state b = 1)
        flux_ab = np.array([expected[2] for expected in fluxes[:,self.state]])
        # CIs in rate (s^-1) format (divided by tau)
        ci_lb_ab = np.array([expected[3] for expected in fluxes[:,self.state]]) * (1/self.tau)
        ci_ub_ab = np.array([expected[4] for expected in fluxes[:,self.state]]) * (1/self.tau)

        # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
        state_pops = np.array(h5["state_pop_evolution"])
        # state A = label 0, state B = label 1
        state_pop_a = np.array([expected[2] for expected in state_pops[:,0]])
        #state_pop_b = np.array([expected[2] for expected in state_pops[:,1]])
        # 2 different approaches here, can norm by state_pop_a (sum of weights in a)
        # but since 2 state system, could also use 1 - state_pop_b since all not in b are in a
        flux_ab = flux_ab / state_pop_a
        #flux_ab = flux_ab / (1 - state_pop_b)

        # convert from tau^-1 to seconds^-1
        rate_ab = flux_ab * (1/self.tau)
        
        return rate_ab, ci_lb_ab, ci_ub_ab

    def plot_rate(self, title=None, moltime=True):
        """
        Plot the rate constant = target flux evolution AB / P_A 

        Returns
        -------
        rate_ab : ndarray
            Array of rates from A -> B in seconds^-1.
        """
        rate_ab, ci_lb_ab, ci_ub_ab = self.extract_rate()

        # WE iterations
        iterations = np.arange(0, len(rate_ab), 1)
        if moltime:
            # convert s to ns
            #self.tau *= 10**-9
            # multiply by tau
            #iterations = np.multiply(iterations, self.tau)
            # convert to ns
            #iterations = np.divide(iterations, 10**9)

            iterations = np.multiply(iterations, 50)
            iterations = np.divide(iterations, 1000)

        if self.units == "mfpts":
            mfpt_ab = 1 / rate_ab
            self.ax.plot(iterations, mfpt_ab, label=self.label)
            #ax.fill_between(iterations, mfpt_ab - (1/ci_lb_ab), mfpt_ab + (1/ci_ub_ab), alpha=0.5)
            self.ax.set_ylabel("MFPT ($s$)")
        elif self.units == "rates":
            self.ax.plot(iterations, rate_ab, label=self.label)
            self.ax.fill_between(iterations, rate_ab - ci_lb_ab, rate_ab + ci_ub_ab, alpha=0.5)
            self.ax.set_ylabel("Rate Constant ($s^{-1}$)")

        if moltime:
            self.ax.set_xlabel(r"Molecular Time (ns)")
        else:
            self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        
        self.ax.set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        self.ax.set_title(title)

        return rate_ab

    def plot_statepop(self):
        """
        Plot the state populations
        """
        # read in direct.h5 file
        h5 = h5py.File(f"{self.scheme}/direct.h5", "r")

        # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
        state_pops = np.array(h5["state_pop_evolution"])
        # state A = label 0, state B = label 1
        state_pop_a = np.array([expected[2] for expected in state_pops[:,0]])
        state_pop_b = np.array([expected[2] for expected in state_pops[:,1]])

        # WE iterations
        iterations = np.arange(0, len(state_pop_a), 1)

        # plot both state population evolutions
        self.ax.plot(iterations, state_pop_a, label="State A")
        self.ax.plot(iterations, state_pop_b, label="State B")
        self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        self.ax.set_ylabel("State Population")

        return state_pop_a, state_pop_b

    def plot_exp_vals(self, ax=None, f_range=False):
        """
        f_range : bool
            Set to True to use mark 25-67 s^-1 as the k_D1D2 rate.
        """
        if ax is None:
            ax = self.ax
        if self.units == "rates":
            if f_range:
                # DTY 19F rates of 25-67 for k_D1D2
                ax.axhspan(25, 67, alpha=0.25, color="grey", label="NMR k$_{D1D2}$")
            else:
                # D1-->D2 ~ 20-50, D2-->D1 ~ 100-150
                ax.axhline(150, color="k", ls="--", label="k$_{D2D1}$")
                ax.axhline(25, color="red", ls="--", label="k$_{D1D2}$")
        elif self.units == "mfpts":
            # converted to mfpt = 1 / rate
            ax.axhline(1/150, color="k", ls="--", label="MFPT$_{D2D1}$")
            ax.axhline(1/25, color="red", ls="--", label="MFPT$_{D1D2}$")
        else:
            raise ValueError(f"You put {self.units} for unit, which must be `mfpts` or `rates`.") 

    # TODO: make this update scheme more automatically
    def plot_multi_def_rates(self, ver="v00", def_subplot=None):
        ### get rates for multiple state definitions (WE c2x 4b)
        final_rates = []
        ogscheme = self.scheme
        for angle in self.states:
            self.scheme = f"{ogscheme}/{self.prefix}{angle}_{ver}"
            self.label = f"> {angle}°"
            rates = self.plot_rate(title=f"2KOD {ver}")
            # add 2 item list: angle | final rate value
            final_rates.append(rates[-1])

        self.plot_exp_vals()

        plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        # TODO: separate into different methods
        # plot the rates at various state definitions
        if def_subplot:
            ax2 = self.fig.add_subplot(def_subplot)
            ax2.plot(self.states, final_rates, color="k")
            ax2.set_xlabel("Angle State Definition", fontsize=11, labelpad=4, color="grey")
            self.plot_exp_vals(ax2)
            ax2.tick_params(colors="grey")
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        #self.ax.set_ylim(10**-9, 10**10)

        self.fig.tight_layout()
        if self.savefig:
            plt.savefig(f"figures/{ogscheme}_{self.prefix}_{ver}multi_state_rates.png", 
                        dpi=300, transparent=True)

    def plot_std_error_rate_reps(self, iterations=500, reps=3, def_subplot=None, title=None):
        """
        ### get rates with std errors for multiple state definitions
        Make a plot of multiple replicates and std err for each tstate def
            Iteration (X) VS Rate (Y) with std error
        And maybe later a plot of multiple state defs (X) vs rate(Y) with error
            Maybe both in one plot again?
        """
        rate_evolutions = np.zeros(shape=(iterations, self.n_states, reps))
        final_rate_scalers = np.zeros(shape=(self.n_states, reps))
        ogscheme = self.scheme
        # for v00, v01, v02, (threshv00)
        for vi, ver in enumerate(range(0,reps)):
        #for vi, ver in enumerate(["v00", "v01", "v02", "threshv00"]):
            for ai, angle in enumerate(self.states):
                self.scheme = f"{ogscheme}/{self.prefix}{angle}_v0{ver}"
                #print(self.scheme)
                self.label = f"> {angle}°"
                # calc and fill out rate
                rate_ab, ci_lb_ab, ci_ub_ab = self.extract_rate()
                # make 3d array of rates: rows = tau, columns = angle, depth = version
                rate_evolutions[:, ai, vi] = rate_ab[:iterations]
                #rate_evolutions[:, ai, vi] = rate_ab
                # array for final rate value
                final_rate_scalers[ai, vi] = rate_ab[-1]

        # plot avg and std error of rate evo
        avg_rate_evo = np.average(rate_evolutions, axis=2)
        stdev_rate_evo = np.std(rate_evolutions, axis=2)
        sterr_rate_evo = stdev_rate_evo / np.sqrt(rate_evolutions.shape[2])
        for ai, angle in enumerate(self.states):
            avg = avg_rate_evo[:,ai]
            err = sterr_rate_evo[:,ai]
            self.ax.plot(avg, label=f"< {angle}°")
            self.ax.fill_between(np.arange(0, iterations, 1), avg - err, avg + err, alpha=0.25)

        self.plot_exp_vals()
        plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        # plt.legend(plot, [f"<{i}°" for i in states],
        #                 loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        # state definition subplot
        if def_subplot:
            # plot the rates at various state definitions
            ax2 = self.fig.add_subplot(def_subplot)

            # plot avg and std error of final rates
            avg_final_scaler = np.average(final_rate_scalers, axis=1)
            stdev_final_scaler = np.std(final_rate_scalers, axis=1)
            #print(final_rate_scalers.shape[1])
            sterr_final_scaler = stdev_final_scaler / np.sqrt(final_rate_scalers.shape[1])

            # TODO: table of final avg +/- stderror values
            rate_table = np.zeros((len(self.states), 3))
            rate_table[:,0] = self.states
            rate_table[:,1] = avg_final_scaler
            rate_table[:,2] = sterr_final_scaler
            #np.savetxt("rate_table.csv", rate_table, delimiter="\t")
            import pandas as pd
            df = pd.DataFrame(rate_table, columns=["State", "Average", "Stderr"])
            print(df)

            ax2.plot(self.states, avg_final_scaler, color="k")
            ax2.fill_between(self.states, 
                             avg_final_scaler - sterr_final_scaler, 
                             avg_final_scaler + sterr_final_scaler, 
                             alpha=0.25, color="k")

            ax2.set_xlabel("Angle State Definition", fontsize=11, labelpad=4, color="grey")
            plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
            self.plot_exp_vals(ax2)
            ax2.tick_params(colors="grey")

        self.ax.set_ylabel("Rate Constant ($s^{-1}$)")
        self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        if title:
            self.ax.set_title(title)

        self.fig.tight_layout()
        if self.savefig:
            plt.savefig(f"figures/{ogscheme}_{self.prefix}_{reps}reps_{iterations}i_avg_std.png", 
                        dpi=300, transparent=True)

        return final_rate_scalers

    # To smooth the rate per tstate def plot
    def smooth(self, y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode="same")
        return y_smooth

    def plot_per_state_def(self, ax):
        final_rate_scalers = self.plot_std_error_rate_reps()

        # plot the rates at various state definitions        
        # plot avg and std error of final rates
        avg_final_scaler = np.average(final_rate_scalers, axis=1)
        stdev_final_scaler = np.std(final_rate_scalers, axis=1)
        sterr_final_scaler = stdev_final_scaler / np.sqrt(final_rate_scalers.shape[1])

        ax.plot(self.states, avg_final_scaler, color="k")
        ax.fill_between(self.states, 
                        avg_final_scaler - sterr_final_scaler, 
                        avg_final_scaler + sterr_final_scaler, 
                        alpha=0.25, color="k")

        #ax.set_xlabel("Angle State Definition", fontsize=16, labelpad=4)
        ax.set_xlabel("Angle State Definition")
        self.ax = ax
        self.plot_exp_vals()
        #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        x = self.states
        y = avg_final_scaler
        # attempts at fitting a curve
        # window size 9, polynomial order 3
        #yhat = savgol_filter(y, 9, 3)
        #yhat = self.smooth(y, 10)

        # technically dx is constant so don't need diff(x)
        dydx = np.diff(y) / np.diff(x)
        ddydx = np.diff(dydx)
        print(ddydx)
        return dydx, ddydx


######################
### 1d oapdt plots ###
######################
fig, ax = plt.subplots()
#k = Kinetics(scheme=f"ANALYSIS2", state=1, ax=ax, tau=5**-10)
#k = Kinetics(scheme=f"WIPA/TEST", state=1, ax=ax, tau=5**-10)
#k.plot_rate()
k = Kinetics(scheme=f"WIPA/FULLY_UNBOUND2", state=0, ax=ax, tau=5**-10, label="Bound to Unbound")
k.plot_rate()
k = Kinetics(scheme=f"WIPA/TO_STUCK", state=0, ax=ax, tau=5**-10, label="Bound to Stuck")
k.plot_rate()
k = Kinetics(scheme=f"WIPA/STUCK_TO_UNBOUND", state=0, ax=ax, tau=5**-10, label="Stuck to Unbound")
k.plot_rate()

ax.axhline(0.05, color="grey", label="Rate from experiment")
#ax.axhline(0.27, color="grey", label="Rate from experiment")
plt.tight_layout()
plt.ylim(1e-54, 300)

plt.legend()
plt.show()
#fig.savefig("eg5_wt_rates2.png", dpi=300, transparent=True)
