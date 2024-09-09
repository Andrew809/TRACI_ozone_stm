import pandas as pd
import numpy as np

def convert_data_to_nparray(data):
    if isinstance(data, pd.Series):
        return data.values
    elif isinstance(data, pd.DataFrame):
        return data.to_numpy()
    else:
        return data  # If it's neither a Series nor a DataFrame, return the original data

# Example usage:

# Create a sample Pandas Series and DataFrame
sample_series = pd.Series([1, 2, 3, 4, 5])
sample_df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})

# Check if variables are Series or DataFrame and convert to NumPy arrays
numpy_series = convert_data_to_nparray(sample_series)
numpy_df = convert_data_to_nparray(sample_df)

series_as_1xn = numpy_series.reshape(1, -1)  #-1 allows np to infer dimension
series_as_nx1 = numpy_series.reshape(-1, 1)

