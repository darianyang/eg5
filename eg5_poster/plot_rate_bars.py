import matplotlib.pyplot as plt
import numpy as np

plt.style.use("~/github/wedap/wedap/styles/default.mplstyle")
plt.figure(figsize=(4, 4))

'''
eg5_wt/WIPA/FULLY_UNBOUND: 8.529091438037076e-49
eg5_wt/WIPA/TO_STUCK: 4.095901503824797e-06
eg5_wt/WIPA/STUCK_TO_UNBOUND: 2.338942153761648e-37
eg5_mon/WIPA/FULLY_UNBOUND: 2.094437321341923e-44
eg5_mon/WIPA/TO_STUCK: 0.021436501066962755
eg5_mon/WIPA/STUCK_TO_UNBOUND: 5.744106263755286e-38
'''

rates = {"WT k$_1$" : 4.095901503824797e-06, 
         "WT k$_2$" : 2.338942153761648e-37,
         "WT k$_{overall}$" : 8.529091438037076e-49,
         "MB k$_1$" : 0.021436501066962755, 
         "MB k$_2$" : 5.744106263755286e-38,
         "MB k$_{overall}$" : 2.094437321341923e-44,
         }


# Extract keys and values
keys = list(rates.keys())
values = list(rates.values())

custom_colors = ['lightsalmon', 'tomato', 'firebrick', 'lightsteelblue', 'cornflowerblue', 'royalblue']

# Create a bar chart
plt.bar(keys, values, color=custom_colors)

# Rotate x-axis labels for better readability (optional)
plt.xticks(rotation=45, ha='right')

plt.axhline(0.05, color="grey", linestyle="--")
plt.axhline(0.007, color="grey", linestyle="--")

# Add labels and title
plt.ylabel('Rate Constant ($s^{-1}$)')
plt.yscale("log")
plt.ylim(1e-54, 300)

# Show the plot
plt.tight_layout()
plt.savefig("rate_bar.png", dpi=600, transparent=True)
plt.show()

