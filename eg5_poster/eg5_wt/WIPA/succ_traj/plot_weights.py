import pickle
import numpy as np
import matplotlib.pyplot as plt

with open('reassigned.pickle', 'rb') as f:
#with open('output.pickle', 'rb') as f:
    data = pickle.load(f)
# (1286, 2, 10) for succ_traj, frames, pathways dataset
print(len(data))
print(len(data[0]))
#print(data[0][0][-1])
#print(data[:][:][-1])
weights = [i[1][-1] for i in data]
#weights = [j[-1] for j in data for i in j]
print(weights)

plt.hist(weights)
plt.show()
