
import pandas as pd
import numpy as np


# FLOW_name_org
# FLOW_casnumber
# LCIAMethod_location
# LCIAMethod_location_name
# CF
# Unit
#
# CF_Uncertainty_Lower
# CF_Uncertainty_Higher
# FLOW_class0
# FLOW_class1
# FLOW_class2
# Species
# LCIAMethod_realm
# LCIAMethod_mathematicalApproach
# Scenario
# CF_derivation











def expand_df(row_series, df, cols_keep_in_df):
    # Create multiindex from the entire row
    #   This will have several entries we drop later
    #   The index repeats the index from the row_series for as many entries are in the df_ffs
    new_multiindex = pd.MultiIndex.from_arrays(np.tile(row_series, (len(df), 1)).transpose(),
                                               names=list(row_series.index))

    # create dataframe from multiindex and the col_values associated with columns to keep
    df_expanded = pd.DataFrame(data=df.loc[:, cols_keep_in_df].to_numpy(),
                               index=new_multiindex,
                               columns=cols_keep_in_df)

    return df_expanded


def prep_df_for_concat(df_target, df_incoming, incoming_cols_keep):
    # df_target and df_incoming share multi-index;
    # df_incoming may have extra (or fewer) multi-level indices.

    # cols_keep_in_df are the (non-multiindex) columns in the df_incoming that we keep

    # check all indices of incoming are present in target:
    if not all(np.in1d(df_target.index.names, df_incoming.index.names)):
        raise KeyError(f'Function <record_df> called with an incoming df with '
                       f'index names not present in target.')

    # fix the index of the incoming...
    # first create a list of those that are not in the target
    mask_matching = np.in1d(df_incoming.index.names, df_target.index.names)
    list_drop = [c for (c, i) in zip(df_incoming.index.names,
                                     mask_matching) if not i]
    # drop unneeded index
    df_incoming.index = df_incoming.index.droplevel(level=list_drop)

    # reorder to match target
    df_incoming = df_incoming.reorder_levels(order=list(df_target.index.names))

    # if there is a set of columns to keep, select those
    df_incoming = df_incoming.loc[:, incoming_cols_keep]

    # add missing columns (concat will do this later,
    #   but we need them now to deal with normalization)
    df_incoming = df_incoming.reindex(columns=df_target.columns)

    return df_incoming