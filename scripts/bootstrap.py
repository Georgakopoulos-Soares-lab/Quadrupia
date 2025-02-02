import numpy as np
import pandas as pd
from tqdm import tqdm

def confidence_interval_bootstrap(df, alpha=0.95, total_resamples = 5000):
    samples = []
    df = df.select_dtypes(include=np.number)
    total_sum = df.sum(axis=0)
    mean = total_sum.mean()
    total_sum = (total_sum / mean).to_frame(name="mean")
    for _ in tqdm(range(total_resamples)):
        resample = df.sample(df.shape[0], replace=True)
        resample_sum = resample.sum(axis=0)
        mean = resample_sum.mean()
        resample_sum = resample_sum / mean
        samples.append(resample_sum)
    samples = pd.DataFrame(samples)
    confidence_intervals = {}
    for col in samples:
        temp = samples[col].sort_values(ascending=True)
        confidence_intervals.update(
                                {col: (np.percentile(temp, 100 * (1-alpha)/2), 
                                           np.percentile(temp, 100 * (alpha+1)/2)
                                        )
                                    }
                                   )
        
    confidence_intervals = pd.DataFrame(confidence_intervals).T
    confidence_intervals.columns = ["inf", "sup"]
    confidence_intervals = pd.concat([total_sum, confidence_intervals], axis=1)
    return confidence_intervals
