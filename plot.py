import numpy as np
import matplotlib.pyplot as plt

def plot_adp(directory):
    prot = np.loadtxt(f"{directory}/rms_protcut_heavy.dat")[::10]
    #complex = np.loadtxt(f"{directory}/rms_complexcut_heavy.dat")
    mon = np.loadtxt(f"{directory}/rms_mon_heavy.dat")[::10]
    adp = np.loadtxt(f"{directory}/rms_adp_heavy.dat")[::10]

    plt.plot(prot[:,0], prot[:,1], label="Protein")
    #plt.plot(complex[:,0], complex[:,1], label="Complex")
    plt.plot(mon[:,0], mon[:,1], label="MON")
    plt.plot(adp[:,0], adp[:,1], label="ADP")
    
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("Eg5 ADP MON")
     
def plot_apo(directory):
    prot = np.loadtxt(f"{directory}/rms_protcut_heavy.dat")[::10]
    #complex = np.loadtxt(f"{directory}/rms_complexcut_heavy.dat")
    #adp = np.loadtxt(f"{directory}/rms_adp_heavy.dat")[::10]
    mon = np.loadtxt(f"{directory}/rms_mon_heavy.dat")[::10]

    plt.plot(prot[:,0], prot[:,1], label="Protein")
    #plt.plot(complex[:,0], complex[:,1], label="Complex")
    #plt.plot(adp[:,0], adp[:,1], label="ADP")
    plt.plot(mon[:,0], mon[:,1], label="MON")
    
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("Eg5 MON No ADP")

plot_adp("prep/eg5_adpmon/md/v00/200ns")
#plot_apo("prep/eg5_adpmon/noadp/md/v00/200ns")
plt.ylim(0,5)
plt.legend()
plt.show()


