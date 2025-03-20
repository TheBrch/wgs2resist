import pandas as pd

from sklearn.model_selection import train_test_split

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import accuracy_score
import joblib

X_bin = pd.read_pickle("data_dense.pkl")
y = pd.read_pickle("labels.pkl")

X_train, X_test, y_train, y_test = train_test_split(X_bin, y, test_size=0.2, random_state=42)

y_train_bin = y_train.values
y_test_bin = y_test.values


models = {
    "gaussian": GaussianProcessClassifier(kernel=RBF(length_scale=1.0)),
    "svm": SVC(C=1.0, kernel="rbf", probability=True),
    "xgboost": joblib.load("xgboost.pkl"),
    "logistic": LogisticRegression(solver="liblinear", penalty="l1")
}

for name, model in models.items():
    print(f"Training {name}...")

    if (name != "xgboost"):
        multi_model = MultiOutputClassifier(model)
        multi_model.fit(X_train, y_train_bin)
        joblib.dump(multi_model, f"{name}.pkl")
    else:
        multi_model = model
    # end if


    y_pred = multi_model.predict(X_test)
    accuracy = accuracy_score(y_test_bin, y_pred)
    print(f"Accuracy of {name}: {accuracy:.4f}")
    
    prob_matrix = np.array([est.predict_proba(X_test)[:, 1] for est in multi_model.estimators_]).T
    print(f"Predicted Probabilities for {name} (rows=samples, columns=labels):")
    print(prob_matrix)
    print("\n")