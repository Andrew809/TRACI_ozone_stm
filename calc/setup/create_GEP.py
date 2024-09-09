

# python match - case is implemented as of v 3.10, so just using switch


from calc.setup import config as cfg
from calc.functions.conversions import convert_data_to_nparray
from calc.inout.read_data import read_column_from_file
import pandas as pd
import numpy as np

# cfg supplies us with choices for processing the GEP:
# str_GEPsetup_normalize_across_receiving = 'normalize across receiving elements'
# str_GEPsetup_normalize_entire_matrix = 'normalize the full [emit x receive] matrix'
# str_GEPsetup_weight_rows_by_vector = 'weight receving rows by vector'

def adjust_GEP_vector(series, adjust_method, num_emit_rows, emit_weights_series=None):
    """

    :param series: a pandas data series (or numpy array) of the GEP values, in the authoritative order (i.e., IDs have been fixed elsewhere).
                    The series must also have had NaNs removed - if not; the normalization will yield NaNs.
    :param adjust_method: a string, set in config, to choose how to make the GEP matrix
    :param num_emit_rows: number of rows in the output matrix.  This should be the number of emission locations in the model.
    :param emit_weights_series: if weighting, a pandas series (or numpy array) of the weights for each emitting row.  Should be in authoritative order.
    :return: a GEP matrix as numpy array
    """

    if not isinstance(series, pd.Series):
        series = series.squeeze()

    if not isinstance(series, pd.Series):
        raise KeyError(f'function <adjust_GEP_vector> was not passed a pd.Series or a pd.Dataframe that could be squeezed into a series')
    else:
        array = convert_data_to_nparray(series)
        vector = array.reshape(1, -1)  #-1 allows numpy to infer the number

        if adjust_method == cfg.str_GEPsetup_normalize_across_receiving:
            norm_vector = vector / np.sum(vector)
            matrix = np.tile(norm_vector, (num_emit_rows, 1))

        elif adjust_method == cfg.str_GEPsetup_normalize_entire_matrix:
            raw_matrix = np.tile(vector, (num_emit_rows, 1))
            matrix = raw_matrix / np.sum(raw_matrix)

        elif adjust_method == cfg.str_GEPsetup_weight_rows_by_vector:
            if emit_weights_series is None:
                raise KeyError(f'function <adjust_GEP_vector> asked to weight rows by vector, but not given the weighting vector')
            else:
                emit_weights = convert_data_to_nparray(emit_weights_series)
                emit_weights = emit_weights.reshape(-1, 1)  # -1 allows numpy to infer the number
                normalized_weights = emit_weights / np.sum(emit_weights)
                #raw_matrix = np.tile(vector, (num_emit_rows, 1))
                # Multiply each row of the data_matrix by the corresponding value in the weights column
                #matrix = raw_matrix * normalized_weights
                matrix = normalized_weights * vector

        return matrix

# testing: 3 emission regions, 5 receiving
test_GEP_values = pd.Series([0.1, 0.2, 0.3, 0.4, 0.5])
test_GEP_weights = pd.Series([50, 20, 10])
test_num_rows = test_GEP_weights.shape[0]
#
# print(f'GEP with {cfg.str_GEPsetup_normalize_across_receiving}')
# test_mat = adjust_GEP_vector(test_GEP_values, cfg.str_GEPsetup_normalize_across_receiving, test_num_rows)
# print(test_mat)
#
# print(f'GEP with {cfg.str_GEPsetup_normalize_entire_matrix}')
# test_mat = adjust_GEP_vector(test_GEP_values, cfg.str_GEPsetup_normalize_entire_matrix, test_num_rows)
# print(test_mat)
#
