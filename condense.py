from sklearn.feature_selection import VarianceThreshold

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split

from xgboost import XGBClassifier
from sklearn.multioutput import MultiOutputClassifier
import joblib

X_bin = pd.read_pickle("snps_bin.pkl")
y = pd.read_pickle("labels.pkl")

X_train, X_test, y_train, y_test = train_test_split(X_bin, y, test_size=0.2, random_state=42)

y_train_bin = y_train.values
y_test_bin = y_test.values

selector = VarianceThreshold(threshold=0.01)
X_thresh = selector.fit_transform(X_bin)

X_thresh = pd.DataFrame(X_thresh)
corr_matrix = X_thresh.corr().abs()
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
to_drop = [column for column in upper.columns if any(upper[column] > 0.9)]

X_reduced = X_thresh.drop(columns=to_drop).values

xgb = XGBClassifier(
    n_estimators=100,
    max_depth=3,
    learning_rate=0.1,
    use_label_encoder=False,
    eval_metric="logloss"
)
xgb_multi = MultiOutputClassifier(xgb)
xgb_multi.fit(X_train, y_train_bin)
joblib.dump(xgb_multi, "xgboost.pkl")

importances = np.mean([est.feature_importances_ for est in xgb_multi.estimators_], axis=0)

threshold = np.mean(importances)
selected_features_mask = importances > threshold
X_final = X_reduced[:, selected_features_mask]

print("Original number of features:", X_bin.shape[1])
print("After variance and correlation filtering:", X_reduced.shape[1])
print("After XGBoost feature selection:", X_final.shape[1])

X_final.to_pickle("data_dense.pkl")