import numpy as np
import pandas as pd

data = {
    "Accuracy": [0.964, 0.978, 0.975, 0.977, 0.971],
    "MCC": [0.813, 0.887, 0.870, 0.883, 0.852],
    "Precision": [0.841, 0.915, 0.905, 0.887, 0.866],
    "Recall": [0.826, 0.884, 0.864, 0.904, 0.871],
    "F1": [0.833, 0.899, 0.884, 0.895, 0.868]
}
df = pd.DataFrame(data)
means = df.mean()
ses = df.std() / np.sqrt(len(df))
summary = pd.DataFrame({
    "Mean": means.round(3),
    "[+-] SE": ses.round(3)
})
print(summary)
