import numpy as np
import matplotlib.pyplot as plt

def plot_timeseries(directory):
    prot = np.loadtxt(f"{directory}/rms_protcut_heavy.dat")[::10]
    #complex = np.loadtxt(f"{directory}/rms_complexcut_heavy.dat")
    adp = np.loadtxt(f"{directory}/rms_adp_heavy.dat")[::10]
    mon = np.loadtxt(f"{directory}/rms_mon_heavy.dat")[::10]

    plt.plot(prot[:,0], prot[:,1], label="Protein")
    #plt.plot(complex[:,0], complex[:,1], label="Complex")
    plt.plot(adp[:,0], adp[:,1], label="ADP")
    plt.plot(mon[:,0], mon[:,1], label="MON")
    
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
     

#plot_timeseries("eg5_adp/md_test2022/v00/equil")
plot_timeseries("eg5_adpmon/md/v00/equil")
#plt.title("2022 Test Eg5 ADP")
plt.title("Eg5 ADP MON")
plt.legend()
plt.show()


