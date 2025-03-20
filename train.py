import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import accuracy_score
import joblib

train_tsv = "training_matrix.tsv"
sus_tsv = "susceptibility.tsv"
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

# REDUCE THE NUMBER OF PARAMETERS

X_train, X_test, y_train, y_test = train_test_split(X_bin, y, test_size=0.2, random_state=42)

lb = LabelBinarizer()
y_train_bin = lb.fit_transform(y_train.values)
y_test_bin = lb.transform(y_test.values)

models = {
    "Gaussian Process": GaussianProcessClassifier(kernel=RBF(length_scale=1.0)),
    "SVM": SVC(C=1.0, kernel="rbf", probability=True),
    # "XGBoost": XGBClassifier(n_estimators=100, max_depth=3, learning_rate=0.1, use_label_encoder=False, eval_metric="logloss"),
    # "Logistic Regression": LogisticRegression(solver="liblinear", penalty="l1")
}

for name, model in models.items():
    print(f"Training {name}...")

    multi_model = MultiOutputClassifier(model)
    multi_model.fit(X_train, y_train_bin)
    
    joblib.dump(multi_model, f"{name}.pkl")

    y_pred = multi_model.predict(X_test)
    accuracy = accuracy_score(y_test_bin, y_pred)
    print(f"Accuracy of {name}: {accuracy:.4f}")
    
    prob_matrix = np.array([est.predict_proba(X_test)[:, 1] for est in multi_model.estimators_]).T
    print(f"Predicted Probabilities for {name} (rows=samples, columns=labels):")
    print(prob_matrix)
    print("\n")