import wedap
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import numpy as np

# load h5 file into pdist class
data = wedap.H5_Pdist("west.h5", data_type="average")

# extract data arrays (can be pcoord or any aux data name)
X = data.get_total_data_array("pcoord", 0)
Y = data.get_total_data_array("pcoord", 1)
Z = data.get_total_data_array("pcoord", 2)

# put X and Y together column wise
XY = np.hstack((X,Y,Z))

# scale data
scaler = StandardScaler()
XY = scaler.fit_transform(XY)

# use 10x less data for easier plotting
#XY = XY[::10,:]

# cluster pdist using weighted k-means
clust = KMeans(n_clusters=18, n_init='auto').fit(XY)

# create plot base
#fig, ax = plt.subplots()
#plot_object = ax.scatter(XY[:,0], XY[:,1], c=clust.labels_, s=1, cmap="tab20")
#plt.colorbar(plot_object)
#plt.show()

# train KNN model to use on other system
from sklearn import neighbors
knn = neighbors.KNeighborsClassifier(n_neighbors=5)
knn.fit(data_train, labels_train)
import pickle
pickle.dump(knn, open("knn.pickle", 'wb'))

# save as new dataset in updated west.h5 file 
wedap.H5_Pdist("west.h5", data_type="average", Xname=clust.labels_, 
               H5save_out="west_labeled.h5", Xsave_name="labels").make_new_h5()
