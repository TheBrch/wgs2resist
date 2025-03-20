import pandas as pd
import numpy as np
import yaml

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

train_tsv = "training_matrix.tsv"
sus_tsv = config["suscept_table"]
sus = pd.read_csv(sus_tsv, sep="\t")
data = pd.read_csv(train_tsv, sep="\t")

antibiotics = [col for col in sus.iloc[:, 0].to_list() if col in data.columns]

y = data[antibiotics]
X = data.drop(columns=antibiotics)

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

X_bin.to_pickle("snps_bin.pkl")
y.to_pickle("labels.pkl")
