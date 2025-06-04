# Copied from gudhi.CoverComplex.fit (version 3.10.1)
# Only change: replaced np.mean with np.nanmean

import numpy as np
import itertools
from gudhi.simplex_tree import SimplexTree
from sklearn.cluster import AgglomerativeClustering
from gudhi.cover_complex import CoverComplex

from sklearn import __version__ as sklearn_version
from sklearn.utils.fixes import parse_version
agglomerative_clustering_metric = parse_version(sklearn_version) >= parse_version('1.2.0')

def fit(self, X, y=None, filters=None, colors=None):
        """
        Fit the MapperComplex class on a point cloud or a distance matrix: compute the Mapper complex and store it in a simplex tree called `simplex_tree_`.

        Parameters
        ----------
            X : numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix
                input point cloud or distance matrix.
            y : n x 1 array
                point labels (unused).
            filters : list of lists or numpy array of shape (num_points) x (num_filters) 
                filter functions (sometimes called lenses) used to compute the cover. Each column of the numpy array defines a scalar function defined on the input points.
            colors : list of lists or numpy array of shape (num_points) x (num_colors)
                functions used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
        """

        if self.resolutions is not None:
            self.resolutions = np.array(self.resolutions)
            if len(self.resolutions.shape) == 0:
                self.resolutions = np.broadcast_to(self.resolutions, 1)
        if self.gains is not None:
            self.gains = np.array(self.gains)
            if len(self.gains.shape) == 0:
                self.gains = np.broadcast_to(self.gains, 1)
        if self.filter_bnds is not None:
            self.filter_bnds = np.array(self.filter_bnds)

        self.filters, self.colors = filters, colors

        if self.filters is None:
            if self.input_type == "point cloud":
                self.filters = X[:,0:1]
            elif self.input_type == "distance matrix":
                self.filters = X.max(axis=0)[:,None]
        else:
            if isinstance(self.filters, np.ndarray) == False:
                self.filters = np.array(self.filters).T

        if self.colors is None:
            if self.input_type == "point cloud":
                self.colors = X[:,0:1]
            elif self.input_type == "distance matrix":
                self.colors = X.max(axis=0)[:,None]
        else:
            if isinstance(self.colors, np.ndarray) == False:
                self.colors = np.array(self.colors).T

        if len(self.filters.shape) == 1: # if self.filters is a 1D filter, convert it to an array of shape [n,1]
            self.filters = np.reshape(self.filters, [len(X),1])
        if len(self.colors.shape) == 1: # if self.colors is a 1D filter, convert it to an array of shape [n,1]
            self.colors = np.reshape(self.colors, [len(X),1])
      
        num_pts, num_filters = self.filters.shape[0], self.filters.shape[1]

        # If some filter limits are unspecified, automatically compute them
        if self.filter_bnds is None:
            self.filter_bnds = np.hstack([np.min(self.filters, axis=0)[:,np.newaxis], np.max(self.filters, axis=0)[:,np.newaxis]])

        # If some resolutions are not specified, automatically compute them
        if self.gains is None:
            self.gains = np.broadcast_to(1/3, num_filters)
        if self.resolutions is None or self.clustering is None:
            delta, resolutions = self.get_optimal_parameters_for_agglomerative_clustering(X=X, beta=self.beta, C=self.C, N=self.N)
            if self.clustering is None:
                self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, **{
        "metric" if agglomerative_clustering_metric else "affinity":
        "euclidean" if self.input_type == "point cloud" else "precomputed"
                                                        })
            if self.resolutions is None:
                self.resolutions = np.multiply(resolutions, 1./self.gains)
                self.resolutions = np.array([int( (self.filter_bnds[ir,1]-self.filter_bnds[ir,0])/r) for ir, r in enumerate(self.resolutions)])

        # Initialize attributes
        self.simplex_tree_, self.node_info_ = SimplexTree(), {}

        if np.all(self.gains < .5):

            # Compute which points fall in which patch or patch intersections
            interval_inds, intersec_inds = np.empty(self.filters.shape), np.empty(self.filters.shape)
            for i in range(num_filters):
                f, r, g = self.filters[:,i], self.resolutions[i], self.gains[i]
                min_f, max_f = self.filter_bnds[i,0], np.nextafter(self.filter_bnds[i,1], np.inf)
                interval_endpoints, l = np.linspace(min_f, max_f, num=r+1, retstep=True)
                intersec_endpoints = []
                for j in range(1, len(interval_endpoints)-1):
                    intersec_endpoints.append(interval_endpoints[j] - g*l / (2 - 2*g))
                    intersec_endpoints.append(interval_endpoints[j] + g*l / (2 - 2*g))
                interval_inds[:,i] = np.digitize(f, interval_endpoints)
                intersec_inds[:,i] = 0.5 * (np.digitize(f, intersec_endpoints) + 1)

            # Build the binned_data map that takes a patch or a patch intersection and outputs the indices of the points contained in it
            binned_data = {}
            for i in range(num_pts):
                list_preimage = []
                for j in range(num_filters):
                    a, b = interval_inds[i,j], intersec_inds[i,j]
                    list_preimage.append([a])
                    if b == a:
                        list_preimage[j].append(a+1)
                    if b == a-1:
                        list_preimage[j].append(a-1)
                list_preimage = list(itertools.product(*list_preimage))
                for pre_idx in list_preimage:
                    try:
                        binned_data[pre_idx].append(i)
                    except KeyError:
                        binned_data[pre_idx] = [i]

        else:

            # Compute interval endpoints for each filter
            l_int, r_int = [], []
            for i in range(num_filters):
                L, R = [], []
                f, r, g = self.filters[:,i], self.resolutions[i], self.gains[i]
                min_f, max_f = self.filter_bnds[i,0], np.nextafter(self.filter_bnds[i,1], np.inf)
                interval_endpoints, l = np.linspace(min_f, max_f, num=r+1, retstep=True)
                for j in range(len(interval_endpoints)-1):
                    L.append(interval_endpoints[j]   - g*l / (2 - 2*g))
                    R.append(interval_endpoints[j+1] + g*l / (2 - 2*g))
                l_int.append(L)
                r_int.append(R)

            # Build the binned_data map that takes a patch or a patch intersection and outputs the indices of the points contained in it
            binned_data = {}
            for i in range(num_pts):
                list_preimage = []
                for j in range(num_filters):
                    fval = self.filters[i,j]
                    start, end = int(min(np.argwhere(np.array(r_int[j]) >= fval))), int(max(np.argwhere(np.array(l_int[j]) <= fval)))
                    list_preimage.append(list(range(start, end+1)))
                list_preimage = list(itertools.product(*list_preimage))
                for pre_idx in list_preimage:
                    try:
                        binned_data[pre_idx].append(i)
                    except KeyError:
                        binned_data[pre_idx] = [i]

        # Initialize the cover map, that takes a point and outputs the clusters to which it belongs
        cover, clus_base = [[] for _ in range(num_pts)], 0

        # For each patch
        for preimage in binned_data:

            # Apply clustering on the corresponding subpopulation
            idxs = np.array(binned_data[preimage])
            if len(idxs) > 1:
                clusters = self.clustering.fit_predict(X[idxs,:]) if self.input_type == "point cloud" else self.clustering.fit_predict(X[idxs,:][:,idxs])
            elif len(idxs) == 1:
                clusters = np.array([0])
            else:
                continue

            # Collect various information on each cluster
            num_clus_pre = np.max(clusters) + 1
            for clus_i in range(num_clus_pre):
                node_name = clus_base + clus_i
                subpopulation = idxs[clusters == clus_i]
                self.node_info_[node_name] = {}
                self.node_info_[node_name]["indices"] = subpopulation
                self.node_info_[node_name]["size"] = len(subpopulation)
                self.node_info_[node_name]["colors"] = np.nanmean(self.colors[subpopulation,:], axis=0)
                self.node_info_[node_name]["patch"] = preimage

            # Update the cover map
            for pt in range(clusters.shape[0]):
                node_name = clus_base + clusters[pt]
                if clusters[pt] != -1 and self.node_info_[node_name]["size"] >= self.min_points_per_node:
                    cover[idxs[pt]].append(node_name)

            clus_base += np.max(clusters) + 1

        # Insert the simplices of the Mapper complex
        for i in range(num_pts):
            self.simplex_tree_.insert(cover[i])

        return self