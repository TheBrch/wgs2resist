from sklearn.feature_selection import VarianceThreshold
import numpy as np
import pandas as pd
import os
import sys
import logging
import gc

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


def get_correlated(X, threshold=0.9, chunksize=1000):
    n_features = X.shape[1]
    to_drop = set()
    hi_corr = {}

    def drop_it(i, j, c):
        if j not in to_drop:
            to_drop.add(j)
            feat_i = X.columns[i]
            feat_j = X.columns[j]
            if feat_j not in hi_corr:
                hi_corr[feat_j] = {}
            hi_corr[feat_j][feat_i] = round(c, 3)

    logging.info("Converting data to NumPy...")
    X_np = X.values

    logging.info("Iterating through features...")

    for i_a in range(0, n_features, chunksize):
        i_b = min(i_a + chunksize, n_features)
        n_i = i_b - i_a

        X_i = X_np[:, i_a:i_b]

        for j_a in range(i_a, n_features, chunksize):
            j_b = min(j_a + chunksize, n_features)
            n_j = j_b - j_a

            X_j = X_np[:, j_a:j_b]

            corr = np.corrcoef(X_i.T, X_j.T)

            if i_a == j_a:
                sub_corr = corr[:n_i, :n_j]
                for sub_i in range(n_i):
                    for sub_j in range(sub_i + 1, n_j):
                        corr_val = sub_corr[sub_i, sub_j]
                        if abs(corr_val) >= threshold:
                            drop_it(i_a + sub_i, j_a + sub_j, corr_val)
            else:
                sub_corr = corr[:n_i, n_i : n_i + n_j]
                for sub_i in range(n_i):
                    up_i = i_a + sub_i
                    if up_i in to_drop:
                        continue
                    for sub_j in range(n_j):
                        corr_val = sub_corr[sub_i, sub_j]
                        if abs(corr_val) >= threshold:
                            up_j = j_a + sub_j
                            drop_it(up_i, up_j, corr_val)
            del corr, sub_corr, X_j
            gc.collect()
        del X_i
        gc.collect()

        if (i_a // chunksize + 1) % 10 == 0:
            logging.info(
                f"Processed {i_b}\t/{n_features} features. Found {len(to_drop)} strongly correlated features."
            )
    to_drop_sorted = sorted(list(to_drop))
    logging.info(f"Dropped number of features: {len(to_drop_sorted)}")

    if hi_corr:
        logging.info(f"Exporting highly correlated features...")
        corr_data = []
        for feature, corrs in hi_corr.items():
            row = {"index": feature}
            row.update(corrs)
            corr_data.append(row)
        corr_df = pd.DataFrame(corr_data)
        corr_df.to_feather(
            os.path.join(
                "results", "condensed_data", f"{antibiotic_name}_hicorr.feather"
            )
        )
        logging.info("Highly correlated features exported.")
    return to_drop_sorted


logging.info("Loading feather file...")
X_bin = pd.read_feather(X_bin_file)

rownames_col = X_bin.pop("rownames")
logging.info(f"Original number of features: {X_bin.shape[1]}")

selector = VarianceThreshold(threshold=0.05)  # try 0.1
X_thresh = pd.DataFrame(
    selector.fit_transform(X_bin), columns=X_bin.columns[selector.get_support()]
)
del X_bin
gc.collect()
logging.info(f"After variance filtering: {X_thresh.shape[1]}")


logging.info(f"Separating highly correlated features...")
to_drop_ind = get_correlated(X_thresh)
logging.info(f"Found {len(to_drop_ind)} highly correlated features")

to_drop_names = X_thresh.columns[to_drop_ind]
X_reduced = X_thresh.drop(columns=to_drop_names)
del X_thresh
gc.collect()
logging.info(f"After correlation filtering: {X_reduced.shape[1]}")

X_reduced["rownames"] = rownames_col

X_reduced.to_pickle(os.path.join("results", "condensed_data", f"{antibiotic_name}.pkl"))
logging.info(f"Exported condensed data. Done.")
