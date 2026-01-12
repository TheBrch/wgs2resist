# import cuml.accel
# cuml.accel.install()

# import cudf.pandas
# cudf.pandas.install()

import numpy as np
import pandas as pd
import yaml
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    auc,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.svm import SVC
from sklearn.base import BaseEstimator, ClassifierMixin
import joblib
import logging
import os
import sys
import xgboost as xgb
from imblearn.over_sampling import SMOTE


class WeightedEnsembleClassifier(BaseEstimator, ClassifierMixin):
    def __init__(self, estimators, fitted_estimators):
        self.estimators = estimators
        self.fitted_estimators = fitted_estimators

    def fit(self, X=None, y=None, groups=None):
        from scipy.optimize import minimize
        from sklearn.utils.validation import check_X_y
        from sklearn.model_selection import cross_val_predict
        from imblearn.pipeline import Pipeline
        from imblearn.over_sampling import SMOTE

        if X is None or y is None:
            raise ValueError("X and y are required for setting weights")
        if groups is None:
            raise ValueError("groups are required for StratifiedGroupKFold")
        X, y = check_X_y(X, y)
        self.classes_ = np.unique(y)

        probas = []
        cv = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=42)
        for name, estimator in self.estimators:
            pipe = Pipeline(
                [("smote", SMOTE(random_state=42)), ("classifier", estimator)]
            )
            oof_proba = cross_val_predict(
                pipe,
                X,
                y,
                groups=groups,
                cv=cv,
                method="predict_proba",
            )[:, 1]
            probas.append(oof_proba)
        probas = np.column_stack(probas)

        def objective(w):
            w = np.maximum(w, 0)
            w /= np.sum(w)
            y_pred = probas @ w
            return np.mean((y - y_pred) ** 2)

        w0 = np.ones(len(self.estimators)) / len(self.estimators)
        bounds = [(0, None)] * len(self.estimators)
        res = minimize(objective, w0, bounds=bounds, method="L-BFGS-B")

        w_opt = np.maximum(res.x, 0)
        w_opt /= w_opt.sum()
        self.weights_ = w_opt
        self.estimators_ = list(self.fitted_estimators)

        return self

    def predict_proba(self, X):
        from sklearn.utils.validation import check_is_fitted, check_array

        check_is_fitted(self)
        X = check_array(X)

        probas = []
        for name, estimator in self.estimators_:
            probas.append(estimator.predict_proba(X)[:, 1])

        return sum(w * p for w, p in zip(self.weights_, probas))

    def predict(self, X):
        proba = self.predict_proba(X)
        return (proba >= 0.5).astype(int)

    def score(self, X, y):
        from sklearn.metrics import accuracy_score

        return accuracy_score(y, self.predict(X))


X_bin_file = sys.argv[1]
antibiotic_name = X_bin_file.split("/")[-1].split(".")[0]

os.makedirs(os.path.join("results", "models", antibiotic_name, "stats"), exist_ok=True)

logging.basicConfig(
    filename=os.path.join("results", "models", antibiotic_name, "training.log"),
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def log_exc(exc_type, exc_value, exc_tb):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))


sys.excepthook = log_exc

with open(os.path.join("config", "config.yaml"), "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

config_models = config["models"].split(" ")
models = [*config_models, "wec"]

X_bin = pd.read_pickle(X_bin_file)
sample_ids = X_bin.pop("rownames")
feature_header = X_bin.columns

patients = sample_ids.str.split("_").str[0]
y = pd.read_pickle(sys.argv[2]).values.ravel()

zero_n_one, counts = np.unique(y, return_counts=True)
splitcount = min(min(counts), 5)

label_stats = pd.DataFrame(
    {
        "Label": zero_n_one,
        "Count": counts,
        "Ratio": [f"{(count / sum(counts)):.2f}" for count in counts],
    }
)

fit_models = []


def define_model(name):
    match name:
        case "xgboost":
            early_stop = xgb.callback.EarlyStopping(
                rounds=10, metric_name="logloss", save_best=True, maximize=False
            )
            return xgb.XGBClassifier(
                n_estimators=1000,
                max_depth=10,
                learning_rate=0.1,
                eval_metric="logloss",
                n_jobs=-1,
                # device="cuda",
                callbacks=[early_stop],
            )
        case "gaussian":
            return GaussianProcessClassifier(
                kernel=RBF(length_scale=1.0), n_restarts_optimizer=10
            )
        case "svm":
            return SVC(C=1.0, kernel="rbf", probability=True)
        case "logistic":
            return LogisticRegression(C=1.0, solver="liblinear", penalty="l1")
        case "wec":
            return WeightedEnsembleClassifier(
                estimators=[(e, define_model(e)) for e in config_models],
                fitted_estimators=fit_models,
            )
        case _:
            logging.error("Unknown model name")
            return


logging.info(f"-----{antibiotic_name}-----\n")
logging.info(
    f"Source susceptibility results:\n{label_stats.to_string(index=False)}\n\n"
)


if splitcount > 1:
    correctness = pd.DataFrame({"sample_id": sample_ids})

    data_collection = pd.DataFrame()
    all_pr = pd.DataFrame()
    all_roc = pd.DataFrame()

    skf = StratifiedGroupKFold(n_splits=splitcount, shuffle=True, random_state=42)
    for fold, (train_index, test_index) in enumerate(skf.split(X_bin, y, patients)):
        train_patients = patients[train_index]
        logging.info(f"\n------Fold {fold}------\n")
        logging.info(f"  Train: index={sample_ids.iloc[train_index].tolist()}")
        logging.info(f"         group={train_patients}")
        logging.info(f"  Test:  index={sample_ids.iloc[test_index].tolist()}")
        logging.info(f"         group={patients[test_index]}")

        X_pre_train, X_test = X_bin.iloc[train_index], X_bin.iloc[test_index]
        y_pre_train, y_test = y[train_index], y[test_index]

        if not np.all(np.isin(zero_n_one, y_pre_train)):
            logging.info(f"Fold {fold} training data not diverse, skipping...")
            continue

        sm = SMOTE(random_state=42)
        X_train, y_train = sm.fit_resample(X_pre_train, y_pre_train)

        for name in models:
            logging.info(f"\nCross-validation of {name} - Fold {fold}...\n")
            model = define_model(name)

            if name == "xgboost":
                model.fit(X_train, y_train, eval_set=[(X_test, y_test)])
            elif name == "wec":
                model.fit(X_pre_train, y_pre_train, train_patients)
            else:
                model.fit(X_train, y_train)
            logging.info(f"Model fitted.")

            fit_models.append((name, model))

            if hasattr(model, "coef_"):
                coef = model.coef_
                featurenames = feature_header
            elif hasattr(model, "feature_importances_"):
                coef = model.feature_importances_
                featurenames = feature_header
            elif hasattr(model, "weights_"):
                coef = model.weights_
                featurenames = config_models
            else:
                coef = np.array([False])

            prob_vector = model.predict_proba(X_test)[:, 1]
            fpr, tpr, roc_thresholds = roc_curve(y_test, prob_vector)
            precision, recall, pr_thresholds = precision_recall_curve(
                y_test, prob_vector
            )
            pr_thresholds = np.append(pr_thresholds, min(prob_vector) - 1e-6)

            f1_scores = 2 * precision * recall / (precision + recall)
            f1_scores = np.nan_to_num(f1_scores)

            best_threshold = pr_thresholds[f1_scores.argmax()]

            roc_df = pd.DataFrame(
                {
                    "fpr": fpr,
                    "tpr": tpr,
                    "threshold": roc_thresholds,
                    "fold": fold,
                    "modelname": name,
                }
            )

            prc_df = pd.DataFrame(
                {
                    "precision": precision,
                    "recall": recall,
                    "threshold": pr_thresholds[: len(precision)],
                    "fold": fold,
                    "modelname": name,
                }
            )

            if np.any(coef != False):
                features = pd.DataFrame(
                    {"feature": featurenames, "value": coef.ravel()}
                )
                features.to_csv(
                    os.path.join(
                        "results",
                        "models",
                        antibiotic_name,
                        "stats",
                        f"{name}_f{fold}_features.tsv",
                    ),
                    sep="\t",
                    index=False,
                )

            # y_pred = model.predict(X_test)
            y_pred = (prob_vector >= best_threshold).astype(int)
            y_correct_pred = (y_test == y_pred).astype(int)
            if name not in correctness.columns:
                correctness[name] = 0
            correctness.iloc[test_index, correctness.columns.get_loc(name)] = (
                y_correct_pred
            )
            cm = confusion_matrix(y_test, y_pred, labels=zero_n_one)

            logging.info(
                f"""Fold {fold} confusion matrix:\n{
                    pd.DataFrame(
                        cm,
                        index=[f"Actual {label}" for label in zero_n_one],
                        columns=[f"Predicted {label}" for label in zero_n_one]
                    )
                }"""
            )

            score = model.score(X_test, y_test)
            logging.info(f"Fold {fold} score: {score}")

            raveled_cm = pd.DataFrame([cm.ravel()], columns=["TN", "FP", "FN", "TP"])

            newdf = pd.DataFrame(
                {
                    "Fold": [fold],
                    "Score": [score],
                    "F1": [f1_score(y_test, y_pred)],
                    "ROC_AUC": [auc(fpr, tpr)],
                    "PR_AUC": [auc(recall, precision)],
                    "Precision@Thresh": [precision_score(y_test, y_pred)],
                    "Recall@Thresh": [recall_score(y_test, y_pred)],
                    "Thresh": [best_threshold],
                    "modelname": [name],
                }
            )

            newdf = pd.concat([newdf, raveled_cm], axis=1)
            data_collection = pd.concat([data_collection, newdf], ignore_index=True)
            all_pr = pd.concat([all_pr, prc_df], ignore_index=True)
            all_roc = pd.concat([all_roc, roc_df], ignore_index=True)

        fit_models = []

    tsv_df = {"crossval_results": data_collection, "roc": all_roc, "pr": all_pr}
    for suffix, dataframe in tsv_df.items():
        for name, group_df in dataframe.groupby("modelname"):
            group_df.drop(columns=["modelname"]).to_csv(
                os.path.join(
                    "results",
                    "models",
                    antibiotic_name,
                    "stats",
                    f"{name}_{suffix}.tsv",
                ),
                sep="\t",
                index=False,
            )

    correctness.to_csv(
        os.path.join("results", "models", antibiotic_name, "stats", f"correctness.tsv"),
        sep="\t",
        index=False,
    )
else:
    logging.info(
        f"{name} - Source data contains too few samples of the least populated class, cross-validation not viable."
    )
