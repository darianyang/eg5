import wedap
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

# load h5 file into pdist class
data = wedap.H5_Pdist("west.h5", data_type="average")

# extract weights
#weights = data.get_all_weights()

# extract data arrays (can be pcoord or any aux data name)
X = data.get_total_data_array("pcoord", 0)
Y = data.get_total_data_array("pcoord", 1)
Z = data.get_total_data_array("pcoord", 2)

#print(X.shape)

# put X and Y together column wise
XY = np.hstack((X,Y,Z))
#print(XY.shape)

# scale data
scaler = StandardScaler()
XY = scaler.fit_transform(XY)

# use 10x less data for easier plotting
#XY = XY[::10,:]

# -ln(W/W(max)) weights
#weights_expanded = -np.log(weights/np.max(weights))[::100]

# cluster pdist using weighted k-means
clust = KMeans(n_clusters=18).fit(XY)#, sample_weight=weights_expanded)

# create plot base
#fig, ax = plt.subplots()

# get color labels
#cmap = np.array(["#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628"])
#colors = [cmap[label] for label in clust.labels_.astype(int)]

#plot_object = ax.scatter(XY[:,0], XY[:,1], c=clust.labels_, s=1, cmap="tab20")

# plot on PCs
#pca = PCA(n_components=2)
#PCs = pca.fit_transform(XY)
#ax.scatter(PCs[:,0], PCs[:,1], c=colors, s=1)

# labels
#ax.set_xlabel("PC1")
#ax.set_ylabel("PC2")
#plt.colorbar(plot_object)
#plt.show()

# save as new dataset in west.h5
wedap.H5_Pdist("west.h5", data_type="average", Xname=clust.labels_, 
               H5save_out="west_labeled.h5", Xsave_name="labels").pdist()
