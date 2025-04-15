import cuml.accel
cuml.accel.install()

import cudf.pandas
cudf.pandas.install()

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
import xgboost as xgb
from sklearn.linear_model import LogisticRegression
import joblib
import os
import sys
import logging

X_bin_file = sys.argv[1]
antibiotic_name = X_bin_file.split("/")[-1].split(".")[0]

os.makedirs(f"models/{antibiotic_name}", exist_ok=True)

logging.basicConfig(
    filename=f"models/{antibiotic_name}/training.log", 
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def log_exc(exc_type, exc_value, exc_tb):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))

sys.excepthook = log_exc

X_bin = pd.read_pickle(X_bin_file)
y = pd.read_pickle(sys.argv[2]).values.ravel()

zero_n_one, counts = np.unique(y, return_counts=True)

label_stats = pd.DataFrame({
    'Label': zero_n_one,
    'Count': counts,
    'Ratio': [f"{(count / sum(counts)):.2f}" for count in counts]
})

def define_xgb():
    early_stop = xgb.callback.EarlyStopping(
        rounds=10,
        metric_name='logloss',
        save_best=True,
        maximize=False
    )

    return xgb.XGBClassifier(
        n_estimators=1000,
        max_depth=10,
        learning_rate=0.1,
        eval_metric="logloss",
        n_jobs=-1,
        device="cuda",
        callbacks=[early_stop]
    )

models = {
    "logistic": LogisticRegression(C=1.0, solver="liblinear", penalty="l1"),
    "gaussian": GaussianProcessClassifier(
        kernel=RBF(length_scale=1.0),
        n_restarts_optimizer=10
    ),
    "svm": SVC(C=1.0, kernel="rbf", probability=True),
    "xgboost": define_xgb()
}

logging.info(f"-----{antibiotic_name}-----\n")
logging.info(f"Source susceptibility results:\n{label_stats.to_string(index=False)}\n\n")

for name, model in models.items():
    logging.info(f"Training {name}...")

    if name != "xgboost":
        model.fit(X_bin, y)
        joblib.dump(model, f"models/{antibiotic_name}/{name}.pkl")
        logging.info(f"{name} model exported.")

    best_xgb_score = -1

    splitcount = min(min(counts), 5)

    logging.info(f"\nCross-validation of {name}...\n")
    if splitcount > 1:
        skf = StratifiedKFold(n_splits=splitcount, shuffle=True, random_state=42)
        data_collection = []
        for fold, (train_index, test_index) in enumerate(skf.split(X_bin, y)):
            logging.info(f"{name} - Fold {fold}...\n")

            X_train, X_test = X_bin.iloc[train_index], X_bin.iloc[test_index]
            y_train, y_test = y[train_index], y[test_index]

            if not np.all(np.isin(zero_n_one, y_train)):
                logging.info(f"Fold {fold} training data not diverse, skipping...")
                continue
            
            if name == "xgboost":
                model.fit(X_train, y_train, eval_set=[(X_test, y_test)])
            else:
                model.fit(X_train, y_train)
            logging.info(f"Model fitted.")

            if hasattr(model, "predict_proba"):
                prob_vector = model.predict_proba(X_test)[:, 1]

            y_pred = model.predict(X_test)
            cm = confusion_matrix(y_test, y_pred, labels=zero_n_one)

            logging.info(
                f'''Fold {fold} confusion matrix:\n{
                    pd.DataFrame(
                        cm,
                        index=[f"Actual {label}" for label in zero_n_one],
                        columns=[f"Predicted {label}" for label in zero_n_one]
                    )
                }'''
            )

            score = model.score(X_test, y_test)
            logging.info(f"Fold {fold} score: {score}")

            newrow = np.concatenate(([score], cm.ravel()))
            if hasattr(model, "predict_proba"):
                newrow = np.concatenate((newrow, prob_vector))
            data_collection.append(newrow)
            
            if name == "xgboost":
                if score > best_xgb_score:
                    best_xgb_fold = fold
                    best_xgb_score = score
                    best_xgb_model = model
                model = define_xgb()

        if name == "xgboost":
            joblib.dump(best_xgb_model, f"models/{antibiotic_name}/{name}.pkl")
            logging.info(f"Best {name} model from fold {best_xgb_fold} exported.")

        df = pd.DataFrame(data_collection)
        column_names = ["Score", "TN", "FP", "FN", "TP"]
        df.columns = column_names + list(df.columns[len(column_names):])
        df.to_csv(f"models/{antibiotic_name}/{name}_crossval_results.tsv", sep='\t', index=True)

    else:
        logging.info(f"{name} - Source data contains too few samples of the least populated class, cross-validation not viable.")
