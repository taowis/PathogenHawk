# interpret_model.py

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def interpret_model(model, feature_names, top_n=10):
    # Get feature importances
    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1]

    # Prepare dataframe
    feat_df = pd.DataFrame({
        'Feature': [feature_names[i] for i in indices],
        'Importance': importances[indices]
    }).head(top_n)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.barh(feat_df['Feature'], feat_df['Importance'], color="steelblue")
    plt.xlabel("Feature Importance")
    plt.gca().invert_yaxis()
    plt.title("Top Feature Importances")
    plt.tight_layout()
    plt.show()
