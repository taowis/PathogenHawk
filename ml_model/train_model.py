import argparse
import yaml
from feature_engineering.build_features import build_features
from ml_model.evaluate_model import evaluate_model
from ml_model.interpret_model import interpret_model
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

def split_data(X, y):
    return train_test_split(X, y, test_size=0.2, random_state=42)

def train_model(X_train, y_train, model_type="xgboost"):
    if model_type == "xgboost":
        model = XGBClassifier(
            use_label_encoder=False,
            objective="binary:logistic",
            eval_metric="logloss",
            base_score=0.5  # must be in (0,1)
        )
        model.fit(X_train, y_train)
        return model
    else:
        raise ValueError(f"Unsupported model: {model_type}")

def main(config_path):
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    X, y = build_features(config)
    X_train, X_test, y_train, y_test = split_data(X, y)
    model = train_model(X_train, y_train, model_type=config["ml_model"])
    evaluate_model(model, X_test, y_test)
    interpret_model(model, X.columns)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to YAML config")
    args = parser.parse_args()

    main(args.config)