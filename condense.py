from sklearn.feature_selection import VarianceThreshold
import numpy as np
import pandas as pd
import os
import sys
import logging

os.makedirs("condensed_data", exist_ok=True)

logging.basicConfig(
    filename=f"condensed_data/{antibiotic_name}.log", 
    level=logging.INFO
)

X_bin_file = sys.argv[1]
antibiotic_name = X_bin_file.split("/")[-1].split(".")[0]


X_bin = pd.read_pickle(X_bin_file)

logging.info("Original number of features:", X_bin.shape[1])

selector = VarianceThreshold(threshold=0.1)
X_thresh = pd.DataFrame(
    selector.fit_transform(X_bin),
    columns=X_bin.columns[selector.get_support()]
)

logging.info("After variance filtering:", X_thresh.shape[1])

corr_matrix = X_thresh.corr().abs()
upper = corr_matrix.where(
    np.triu(
        np.ones(corr_matrix.shape),
        k=1
    ).astype(bool)
)

to_drop = [column for column in upper.columns if any(upper[column] > 0.9)]

X_reduced = X_thresh.drop(columns=to_drop)

logging.info("After correlation filtering:", X_reduced.shape[1])

X_reduced.to_pickle(f"condensed_data/{antibiotic_name}.pkl")