import pandas as pd
import numpy as np
import yaml
import os
import sys

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

data_file = sys.argv[1]
antibiotic_name = data_file.split("/")[-1].split(".")[0]
data = pd.read_csv(data_file, sep="\t")

y = data[["Susceptible"]]
X = data.drop(columns=["Susceptible"])

nuc_encoding = {
    'A': [0, 0],
    'G': [0, 1],
    'T': [1, 0],
    'C': [1, 1]
}

def encode_nuc(s):
    return nuc_encoding.get(s, [None, None])

new_cols = []
for col in X.columns:
    encoded = X[col].apply(encode_nuc)
    new_cols.append(X[col].apply(lambda x: encode_nuc(x)[0]).rename(f'{col}_bit1'))
    new_cols.append(X[col].apply(lambda x: encode_nuc(x)[1]).rename(f'{col}_bit2'))

X_bin = pd.concat(new_cols, axis=1)

os.makedirs("binarized_data", exist_ok=True)

X_bin.to_pickle(f"binarized_data/{antibiotic_name}.pkl")
y.to_pickle(f"binarized_data/{antibiotic_name}_lab.pkl")
