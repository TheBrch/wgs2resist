from sklearn.feature_selection import VarianceThreshold
from numba import njit, prange
import numpy as np
import pandas as pd
import os
import sys
import logging

os.makedirs(os.path.join("results", "condensed_data"), exist_ok=True)

X_bin_file = sys.argv[1]
antibiotic_name = X_bin_file.split("/")[-1].split(".")[0]

logging.basicConfig(
    filename=os.path.join("results", "condensed_data", f"{antibiotic_name}.log"),
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def log_exc(exc_type, exc_value, exc_tb):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))


sys.excepthook = log_exc

X_bin = pd.read_feather(X_bin_file)
rownames_col = X_bin.pop("rownames")

logging.info(f"Original number of features: {X_bin.shape[1]}")

selector = VarianceThreshold(threshold=0.05)  # try 0.1
X_thresh = pd.DataFrame(
    selector.fit_transform(X_bin), columns=X_bin.columns[selector.get_support()]
)

logging.info(f"After variance filtering: {X_thresh.shape[1]}")


@njit(parallel=True)
def get_correlated(corr_matrix):
    threshold = 0.9
    n = corr_matrix.shape[0]
    to_drop = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        for j in range(i + 1, n):
            if abs(corr_matrix[i, j]) > threshold:
                to_drop[j] = True
    return np.nonzero(to_drop)[0]


corr_matrix = X_thresh.corr().astype(np.float32)
logging.info(f"Generated correlation matrix.")

to_drop_ind = get_correlated(corr_matrix.values)
logging.info(f"Found {len(to_drop_ind)} highly correlated SNP bits")

to_drop = corr_matrix.columns[to_drop_ind]
hicorr = corr_matrix[to_drop]
hicorr.reset_index(inplace=True)
hicorr.to_feather(
    os.path.join("results", "condensed_data", f"{antibiotic_name}_hicorr.feather")
)
logging.info(f"Exported highly correlated SNP bit info.")

X_reduced = X_thresh.drop(columns=to_drop)
logging.info(f"After correlation filtering: {X_reduced.shape[1]}")

X_reduced["rownames"] = rownames_col

X_reduced.to_pickle(os.path.join("results", "condensed_data", f"{antibiotic_name}.pkl"))
logging.info(f"Exported condensed data. Done.")
