import sklearn.datasets
import sklearn.metrics
from sklearn.model_selection import train_test_split
import xgboost as xgb


def train_breast_cancer(config):
    # Load dataset
    data, labels = sklearn.datasets.load_breast_cancer(return_X_y=True)
    # Split into train and test set
    train_x, test_x, train_y, test_y = train_test_split(
        data, labels, test_size=0.25)
    # Build input matrices for XGBoost
    train_set = xgb.DMatrix(train_x, label=train_y)
    test_set = xgb.DMatrix(test_x, label=test_y)
    # Train the classifier
    results = {}
    print(config)
    bst = xgb.train(
        config,
        train_set,
        evals=[(test_set, "eval")],
        evals_result=results,
        verbose_eval=False)
    return results


if __name__ == "__main__":
    results = train_breast_cancer({
        "objective": "binary:logistic",
        "eval_metric": ["logloss", "error"]
    })
    accuracy = 1. - results["eval"]["error"][-1]
    print(f"Accuracy: {accuracy:.4f}")