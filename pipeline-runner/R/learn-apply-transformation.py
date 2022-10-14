from scanorama import transform
from sklearn.neighbors import NearestNeighbors

# Apply integrated coordinates back to full data.
def apply_transf(unsketched_data, sketched_data, sketched_integrated_data):
    for i, (X_unsketched, X_sketched) in enumerate(zip(unsketched_data, sketched_data)):
        neigh = NearestNeighbors(n_neighbors=3).fit(X_unsketched)
        _, neigh_idx = neigh.kneighbors(X_sketched)

        ds_idxs, ref_idxs = [], []
        for ref_idx in range(neigh_idx.shape[0]):
            for k_idx in range(neigh_idx.shape[1]):
                ds_idxs.append(neigh_idx[ref_idx, k_idx])
                ref_idxs.append(ref_idx)

        X_int = sketched_integrated_data[i]
        bias = transform(X_unsketched, X_int, ds_idxs, ref_idxs, 15, batch_size=1000)

        sketched_integrated_data[i] = X_unsketched + bias

    return sketched_integrated_data