import cuml.accel
cuml.accel.install()

import cudf.pandas
cudf.pandas.install()

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import joblib
import os
import sys
import logging

X_bin_file = sys.argv[1]
antibiotic_name = X_bin_file.split("/")[-1].split(".")[0]

os.makedirs(f"models/{antibiotic_name}", exist_ok=True)

logging.basicConfig(
    filename=f"models/{antibiotic_name}/training.log", 
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def log_exc(exc_type, exc_value, exc_tb):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))

sys.excepthook = log_exc

X_bin = pd.read_pickle(X_bin_file)
y = pd.read_pickle(sys.argv[2]).values

rnd_s = 41
y_train = np.array([0])

while np.all(y_train == y_train[0]):
    rnd_s+=1
    X_train, X_test, y_train, y_test = train_test_split(X_bin, y, test_size=0.2, random_state=rnd_s)
    logging.warning(f"Test split with seed {rnd_s-1} produced an identical label set. Retrying with seed {rnd_s} ...")


models = {
    "logistic": LogisticRegression(C=1.0, solver="liblinear", penalty="l1"),
    "gaussian": GaussianProcessClassifier(kernel=RBF(length_scale=1.0)),
    "svm": SVC(C=1.0, kernel="rbf", probability=True),
    "xgboost": XGBClassifier(
        n_estimators=100,
        max_depth=3,
        learning_rate=0.1,
        use_label_encoder=False,
        eval_metric="logloss",
        n_jobs=-1,
        device="cuda"
    )
}

for name, model in models.items():
    logging.info(f"Training {name}...")

    model.fit(X_train, y_train)

    joblib.dump(model, f"models/{antibiotic_name}/{name}.pkl")

    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    logging.info(f"Accuracy of {name}: {accuracy:.4f}")
    
    if hasattr(model, "predict_proba"):
        prob_vector = model.predict_proba(X_test)[:, 1]
        logging.info(f"Predicted Probabilities for {name}:")
        logging.info(f"{prob_vector}")
    logging.info("\n")