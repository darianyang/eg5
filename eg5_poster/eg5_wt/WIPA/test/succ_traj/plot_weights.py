import pickle
import numpy as np
import matplotlib.pyplot as plt

#with open('reassigned.pickle', 'rb') as f:
with open('output.pickle', 'rb') as f:
    data = pickle.load(f)
# (1286, 2, 10) for succ_traj, frames, pathways dataset
print(len(data))
print(data[0])
#print(data[:,:,-1])


#plt.plot(data[:,-1])
#plt.show()
