import pandas as pd
import numpy as np
import os
import sys

data_file = sys.argv[1]
antibiotic_name = data_file.split("/")[-1].split(".")[0]
data = pd.read_csv(data_file, sep="\t")

rownames = data[["rownames"]]
y = data[["Susceptible"]]
X = data.drop(columns=["Susceptible", "rownames"])

nuc_encoding = {
    "A": [1, 0, 0, 0],
    "G": [0, 1, 0, 0],
    "T": [0, 0, 1, 0],
    "C": [0, 0, 0, 1],
    "-": [0, 0, 0, 0],
}


def encode_nuc(s):
    return nuc_encoding.get(s, [None, None, None, None])


new_cols = [rownames]
for col in X.columns:
    encoded = X[col].apply(encode_nuc)
    new_cols.append(X[col].apply(lambda x: encode_nuc(x)[0]).rename(f"{col}_A"))
    new_cols.append(X[col].apply(lambda x: encode_nuc(x)[1]).rename(f"{col}_G"))
    new_cols.append(X[col].apply(lambda x: encode_nuc(x)[2]).rename(f"{col}_T"))
    new_cols.append(X[col].apply(lambda x: encode_nuc(x)[3]).rename(f"{col}_C"))

X_bin = pd.concat(new_cols, axis=1)

os.makedirs(os.path.join("results", "binarized_data"), exist_ok=True)

X_bin.to_feather(
    os.path.join("results", "binarized_data", f"{antibiotic_name}.feather")
)
y.to_pickle(os.path.join("results", "binarized_data", f"{antibiotic_name}_lab.pkl"))
