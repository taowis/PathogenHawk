import pandas as pd
from sklearn.preprocessing import LabelEncoder

def build_features(config):
    # Load phenotype/label data
    metadata = pd.read_csv(config["phenotype_file"])

    # Load real feature matrix (from config["feature_file"])
    feature_data = pd.read_csv(config["feature_file"], sep="\t")

    # Merge on sample_id
    merged = pd.merge(metadata, feature_data, on="sample_id")

    # Extract all columns except 'sample_id'， label and fluconazole_MIC
    exclude_cols = {"sample_id", config["labels"], "fluconazole_MIC"}
    feature_cols = [col for col in merged.columns if col not in exclude_cols]

    # X = features, y = encoded labels
    X = merged[feature_cols]
    y = LabelEncoder().fit_transform(merged[config["labels"]])

    return X, y