from sklearn.feature_selection import VarianceThreshold
from numba import jit, prange
import numpy as np
import pandas as pd
import os
import sys
import logging

os.makedirs("condensed_data", exist_ok=True)

X_bin_file = sys.argv[1]
antibiotic_name = X_bin_file.split("/")[-1].split(".")[0]

logging.basicConfig(
    filename=f"condensed_data/{antibiotic_name}.log", 
    level=logging.DEBUG
)
def log_exc(exc_type, exc_value, exc_tb):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))

sys.excepthook = log_exc

X_bin = pd.read_pickle(X_bin_file)

logging.info(f"Original number of features: {X_bin.shape[1]}")

selector = VarianceThreshold(threshold=0.1)
X_thresh = pd.DataFrame(
    selector.fit_transform(X_bin),
    columns=X_bin.columns[selector.get_support()]
)

logging.info(f"After variance filtering: {X_thresh.shape[1]}")

@jit(nopython=True, parallel=True)
def get_correlated(corr_matrix, threshold=0.9):
    to_drop = set()
    n = corr_matrix.shape[0]
    for i in prange(n):  # Parallel loop
        for j in range(i + 1, n):
            if corr_matrix[i, j] > threshold:
                to_drop.add(j)
    return sorted(to_drop)

corr_matrix = X_thresh.corr().abs().astype(np.float32)
corr_matrix.to_csv(f"condensed_data/{antibiotic_name}_corr.tsv", sep="\t", index=True)
to_drop = get_correlated(corr_matrix.values)

X_reduced = X_thresh.drop(columns=to_drop)

logging.info(f"After correlation filtering: {X_reduced.shape[1]}")

X_reduced.to_pickle(f"condensed_data/{antibiotic_name}.pkl")