import yaml
from genome_processing import variant_calling
from feature_engineering import build_features
from ml_model import train_model, evaluate_model, interpret_model

def run_pipeline(config_file):
    with open(config_file) as f:
        cfg = yaml.safe_load(f)
    variant_calling.run(cfg)
    X, y = build_features.run(cfg)
    model = train_model.run(X, y, cfg)
    evaluate_model.run(model, X, y)
    interpret_model.run(model, X)

if __name__ == "__main__":
    import sys
    run_pipeline(sys.argv[1])