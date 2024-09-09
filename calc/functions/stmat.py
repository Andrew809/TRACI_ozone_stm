from calc.setup import config as cfg

import pandas as pd
import numpy as np

from calc.functions.cleanup import remove_string_from_array



# global parameters to be used by calling function
s_stmShape_origXnew = 'orig x new'
s_stmShape_newXorig = 'new x orig'
s_stmNormalize_rows = 'row sums = 1' # sum along a row
s_stmNormalize_cols = 'col sums = 1' # sum along a column

const_frac_for_NaN_weight = 1e-6


def create_stm_from_df(dframe, stm_shape, stm_normalize, col_orig_IDs, col_new_IDs, col_isect_area,
                       col_weight_IDs='', col_weight_vals ='', weight_type='', use_area_weight_if_needed=False):
    r"""
    Create a spatial transform matrix, with or without weighting
    An STM has rows and columns corresponding to orig and new spatial geometries.
    These can be in either (orig x new) or (new x orig) configuration, depending on whether the transform is on the receiving or emitting dimension of, e.g., an FF matrix
    Note that  this works from the IDs provided... so if some are missing, the list of IDs (the row and column indices) will be missing these as well.
        One of the calling functions should check this / reindex to have the authoritative list of IDs

    For GLAM acidification, the two 'layouts' we use are
        to transform receiving FF units to, e.g., ecoregions: : 'orig x new', 'row sums = 1', no weighting
        FF' = FF x STM
        to transform emitting FF units to, e.g., countries: 'new x orig', 'col sums = 1', with weighting
        FF' = STM x FF

    #TODO:
    # if weight values are missing (not 0, but missing) for an entire target area, should default to area-weighting

    See example for checking here:
        W:\work\tools\software\excel\adh_samples\weighting for aggregation.xlsx

        Isect_ID	Targ_ID	Source_ID	Weight_ID	Isect_Area	Weight_Val
        1	A	1	alpha	0.4	1
        2	C	1	alpha	0.1	1
        3	C	2	alpha	0.5	1
        4	A	1	beta	0.2	5
        5	C	1	beta	0.2	5
        6	C	2	beta	0.5	5
        7	B	1	beta	0.1	5
        8	B	3	beta	0.4	5
        9	C	3	beta	0.1	5
        10	C	4	beta	0.5	5
        11	B	3	gamma	0.4	10
        12	C	3	gamma	0.1	10
        13	C	4	gamma	0.5	10


    :param dframe: dataframe of intersection geometry output.  Each row is a unique intersection.
    :param stm_shape: see variables declared at top of this module
    :param stm_normalize: see variables declared at top of this module
    :param col_orig_IDs: dataframe column name corresponding to the IDs (names or numbers) of the original spatial units
    :param col_new_IDs:
    :param col_isect_area: dataframe column name with the area of the intersection units
    :param col_weight_IDs: if used...
    :param col_weight_vals: if used, the values (intensive or extensive) used for weighting the transform
    :param weight_type: extensive or intensive
    :param use_area_weight_if_needed: boolean; if false, missing weights (i.e., NaN) are set to very small values.
            This means that a target with no weights uses area-weighting.  See the test excel file noted above.
    :return: a dataframe with index and column labels, and fractions in the cells i,j.
    """

    # for area-weighting as a backup for missing weights:
    # this factor has looked adequate in testing.  Difficult to check in all cases, but easy to confirm by tweaking and comparing results.
    c_factor_for_missing_weights = 1e6


    # TODO - drop all the unnecessary columns?
    # define some names of new columns we calculate
    calc_area_orig = 'Area Orig'
    calc_area_new = 'Area New'
    calc_area_weight = 'Area Weight'
    calc_use_weight_value = 'Use Weight Value'
    calc_areaXweight = 'Area x Weight Value'

    if cfg.bln_debug: print(f'\tRunning function <create_stm_from_df>...')

    df = dframe.copy()

    # we are weighting if a weight type, IDs, and values are provided
    bln_weight = (weight_type != '') & (col_weight_IDs != '') & (col_weight_vals != '')

    # clean up some of the indices
    # we don't need to sort this yet; (unique doesn't sort), and we do a sort at the end
    # returns numpy nd array
    ids_orig = df.loc[:, col_orig_IDs].dropna().unique()
    ids_orig = np.sort(ids_orig)
    ids_new = df.loc[:, col_new_IDs].dropna().unique()
    ids_new = np.sort(ids_new)
    if bln_weight:
        ids_weight = df.loc[:, col_weight_IDs].dropna().unique()
        ids_weight = np.sort(ids_weight)
    # in case there is an empty string in the ids:
    ids_orig = remove_string_from_array(nparray=ids_orig, str_to_remove='')
    ids_new = remove_string_from_array(nparray=ids_new, str_to_remove='')
    if cfg.bln_debug2: print(f'\t\t\t In function <create_stm_from_df>,  unique ids_new, after cleaning (first 20): {ids_new[0:20]}')
    if bln_weight:
        ids_weight = remove_string_from_array(nparray=ids_weight, str_to_remove='')


    # aggfunction = 'sum'  # in other versions of this, could use mean if forcing values to 1, because all weight values set to 1

    # get the total area of orig, new, weight and add to the dataframe:
    # this adds new columns, where each intersection unit has the total area of the orig, new, or weight
    df[calc_area_orig] = df[col_isect_area].groupby(df[col_orig_IDs]).transform('sum')
    df[calc_area_new] = df[col_isect_area].groupby(df[col_new_IDs]).transform('sum')
    if bln_weight:
        df[calc_area_weight] = df[col_isect_area].groupby(df[col_weight_IDs]).transform('sum')
        # also adjust weight for intensive (no adjustment) or extensive (where weight value is divided by total area of the weight)
        if weight_type == cfg.str_weightTypeinten:
            df[calc_use_weight_value] = df[col_weight_vals]
        elif weight_type == cfg.str_weightTypeexten:
            df[calc_use_weight_value] = df[col_weight_vals] / df[calc_area_weight]

        if use_area_weight_if_needed:
            # replace NaNs in the 'Use Weight Value' column with a small number; e.g., 1e-6 of the column min
            # This is funcitonally equivalent to using area-weighting when weights are missing for an entire region,
            # because in the areas with no value weights, all the weights are now equal.
            # In areas that have some valid and some non-valid,
            #   the non-valid weights are replaced with negligible, non-zero values, so they don't affect
            col_min = min(df[calc_use_weight_value])
            replace_weight = col_min * const_frac_for_NaN_weight
            df.loc[:, calc_use_weight_value] = df.loc[:, calc_use_weight_value].fillna(replace_weight)
        else:
            # do nothing; the NaNs stay and could mean that a target with no valid weight values (only NaNs) doesn't have any final values
            pass

        # for each row (which is an intersection), calculate the intersect area * weight value
        df[calc_areaXweight] = df[col_isect_area] * df[calc_use_weight_value]

    else:
        df[calc_areaXweight] = df[col_isect_area]
    #note that the df will have NaNs where the orig and new shapefiles do not overlap.

    # pivot and sum, setting the rows and columns appropriately
    if stm_shape == s_stmShape_origXnew:
        stm = pd.pivot_table(df, index=col_orig_IDs, columns=col_new_IDs, values=calc_areaXweight, aggfunc='sum')
        # pd can introduce NaNs where there is no data at the new/orig intersection.  We fix later.
        # fill it out... even with 'dropna=False' the pivot_table drops rows with na
        if cfg.bln_debug2: print(f'stm after pivot ...\n {stm.head()}')
        stm = stm.reindex(index = ids_orig, columns = ids_new)
        stm.index.set_names(names=col_orig_IDs, level=None, inplace=True)
        stm.columns.set_names(names=col_new_IDs, level=None, inplace=True)

    else:
    #elif stm_shape == s_stmShape_newXorig:
        stm = pd.pivot_table(df, index=col_new_IDs, columns = col_orig_IDs, values = calc_areaXweight, aggfunc ='sum')
        # fill it out... even with 'dropna=False' the pivot_table drops rows with na
        if cfg.bln_debug2: print(f'stm after pivot ...\n {stm.head()}')
        stm = stm.reindex(index = ids_new, columns = ids_orig)
        stm.index.set_names(names=col_new_IDs, level=None, inplace=True)
        stm.columns.set_names(names=col_orig_IDs, level=None, inplace=True)


    if cfg.bln_debug2: print(f'stm after pivot and reindex...\n {stm.head()}')

    # convert to numeric...
    # we have now done this in calling functions, but we repeat here
    # sometimes data were converting to float
    try:
        stm.index = stm.index.astype('int')
        stm.columns = stm.columns.astype('int')
    except (TypeError, ValueError):
        # in testing, got Type and Value errors
        # if we cannot convert to int (e.g., index is strings)
        pass

    # normalize by rows or columns:
    if stm_normalize == s_stmNormalize_rows:
        # we want the row sums (summing horizontally, along the i values) to be 1
        # so sum axis 1 gives a [row x 1] vector, then divide by that along the rows
        stm_normalize = stm.div(stm.sum(axis=1), axis=0)
    else:
    #elif stm_normalize == s_stmNormalize_cols:
        # we want the column sums (summing vertically, along the j values) to be 1
        # so sum axis 0 gives a 1 x column vector, then divide by that
        stm_normalize = stm.div(stm.sum(axis=0),axis=1)

    if cfg.bln_debug2: print(f'stm after normalize ...\n {stm_normalize.head()}')
    if cfg.bln_debug: print(f'\t\t .... done with function <create_stm_from_df>')

    # Fill NaNs: note that there may be NaNs from
    #   a) the initial pivot table (no values at the intersection of the orig and new)
    #   b) the normalization (if the row or col sum is zero, and we divided)
    stm_normalize.fillna(value=0, inplace=True)
    if cfg.bln_debug2: print(f'stm after filling na: number of nas = {stm_normalize.isnull().sum().sum()} ...\n {stm_normalize.head()}')


    return stm_normalize


#
# def create_sdm_from_df2(dframe, col_i, col_j, col_vals='',
#                         col_i_area='', col_j_area='',
#                         col_intersect_area='', divide_area_by='',
#                         force1to1=False, i_values_intensive=True,
#                         normalize_by_axis=''):
#     """
#     Create a spatial distribution matrix, used to spatially translated between
#     two different geometries.  Rows(i) are 'input' and Columns(j) are output.
#     :param dframe: dataframe of intersections; each row has i and j ids,
#     and may contain the area of their intersection
#     :param col_orig_IDs: name of the column with i ids
#     :param col_j:
#     :param col_vals: name of the column with col_values (if not using area fractions)
#     :param col_i_area: name of the column with i areas
#     :param col_j_area:
#     :param col_intersect_area: name of the column with intersection areas
#     :param divide_area_by: do we normalize by i or j (see below for discussion of intensive/extensive)
#     :param force1to1: boolean; if true, put a 1 at each intersection
#     :param i_values_intensive: boolean: are the i col_values intensive?
#     :param normalize_by_axis: force the max col_values per row or column to be 1
#     :return: a dataframe as described above.
#     """
#
#     print(f'\tRunning function <create_sdm_from_df2>...')
#
#     #    % SDMs.  This is an SDM where each row i, col j = is the fraction of air cell i is in lme j
#     #    %     The sum along a row here are <= 1
#     # need to think about what i,j means.
#
#     # TODO: for aggfunction if col_vals are provided,
#     #  could check (drop_duplicates) there are not multiple rows with i,j the same.
#     #  We will have multiple rows in the case of a 3-shapefile intersect,
#     #  but with just 2, this shouldn't be an issue
#
#     df = dframe.copy()
#     bln_proceed = False
#     aggfunction = ''
#
#     # note:
#     # for air, this returns strings: array(['295', '296', '297', ..., '12300', '12301', '12302'], dtype=object)
#
#     # we don't need to sort this yet; (unique doesn't sort) we do a sort at the end
#     # returns numpy nd array
#     ids_i = df.loc[:, col_i].dropna().unique()
#     ids_j = df.loc[:, col_j].dropna().unique()
#
#     ids_i = remove_string_from_array(nparray=ids_i, str_to_remove='')
#     ids_j = remove_string_from_array(nparray=ids_j, str_to_remove='')
#     # print(f'unique ids_j: {ids_j}')
#
#     if col_vals != '':
#         # we can proceed directly; ignore the col_i_areas and col_j_areas
#         bln_proceed = True
#         # we ASSUME that there is only one value per intersection
#         aggfunction = 'sum'
#     else:
#         # we get into an area division (messy) situation:
#         if col_vals == '':
#             # we have to get col_vals through areas or setting to 1
#             col_vals = 'ThisFunctionColVals'  # hard-coded; used internally in this function
#
#             if force1to1:
#                 # will set all intersections to 1
#                 bln_proceed = True
#                 df[col_vals] = 1
#                 # note that when we do the pivot, we'll take the mean
#                 # in case there are multiple intersections, we just want a value of one.
#                 aggfunction = 'mean'
#
#             elif col_intersect_area != '' and (col_i_area != '' or col_j_area != ''):
#                 # elif intersect_area is not blank, and either i_area or j_area not blank
#                 # Calculate area fractions first
#                 bln_proceed = True
#
#                 # # drop the na col_values, which should also remove empty indices
#                 # # i.e., na would occur in an area if
#                 # # there was not an intersection with the i or j in question
#                 # if col_i_area != '':
#                 #     df_ffs.dropna(subset=[col_i_area], inplace=True)
#                 # if col_j_area != '':
#                 #     df_ffs.dropna(subset=[col_j_area], inplace=True)
#
#                 # calculate the col_values as the fraction of j or i in the intersect
#                 # We will use SDM to calculate an overall value for polygon j,
#                 # so this is this is area-based weighting with intensive col_values.
#                 # (i.e., the col_values from col i are not adjusted based on area fractions)
#                 # maybe need to watch out for division by zero
#
#                 if divide_area_by == 'i':
#                     df[col_vals] = df[col_intersect_area] / df[col_i_area]
#                     if not i_values_intensive:
#                         # cannot do area-weighting by i and also have it extensive... I think
#                         bln_proceed = False
#
#                 elif divide_area_by == 'j':
#                     df[col_vals] = df[col_intersect_area] / df[col_j_area]
#
#                     if not i_values_intensive:
#                         # we further reduce the col_vals by the i area fractions
#                         df[col_vals] = (df[col_vals] *
#                                         (df[col_intersect_area] / df[col_i_area]))
#                 else:
#                     bln_proceed = False
#
#                 aggfunction = 'sum'
#
#             else:
#                 # cannot proceed
#                 bln_proceed = False
#                 return 0
#
#     if not bln_proceed:
#         return 0
#     else:
#         sdm = pd.pivot_table(df, index=col_i, columns=col_j, values=col_vals,
#                              aggfunc=aggfunction)
#
#         # fill it out... even with 'dropna=False' the pivot_table drops rows with na
#         sdm = sdm.reindex(index = ids_i, columns = ids_j)
#
#         # convert to numeric...
#         # we have now done this in calling functions, but we repeat here
#         # sometimes data were converting to float
#         try:
#             sdm.index = sdm.index.astype('int')
#             sdm.columns = sdm.columns.astype('int')
#         except (TypeError, ValueError):
#             # in testing, got Type and Value errors
#             # if we cannot convert to int (e.g., index is strings)
#             pass
#
#
#         # sort 'em
#         sdm.sort_index(axis='index', inplace=True)
#         sdm.sort_index(axis='columns', inplace=True)
#
#         sdm.fillna(value=0, inplace=True)
#         # print(sdm.head())
#
#         sdm.index.set_names(names=col_i, level=None, inplace=True)
#         sdm.columns.set_names(names=col_j, level=None, inplace=True)
#
#         if normalize_by_axis == 0:
#             # make it so row sum is 1
#             # (i.e., the sum will add down a column to create a 1 x n vector of 1s)
#             sdm = sdm.div(sdm.sum(axis=0), axis=1)
#
#         elif normalize_by_axis == 1:
#             sdm = sdm.div(sdm.sum(axis=1), axis=0)
#
#             # TODO from matlab version, could check that sums are < 1:
#             # if bln_normalize
#             #     % divide by columns of the A value sums
#             #     m_ab = m_ab ./ valuesA ;
#             #
#             #     % check sums of rows
#             #     rowSums = sum(m_ab,2) ;
#             #
#             #     if ~exist('tol','var')
#             #         tol = 1e-6 ;
#             #     end
#             #
#             #     if any(rowSums > 1 + tol)
#             #         fprintf ('whoops; rows sums are greater than 1 + tol.  Max = %0.6f.\n', max(rowSums)) ;
#             #         bln_error = true ;
#             #     end
#             # end
#         print(f'\t\t... done with function <create_sdm_from_df2>.')
#         return sdm


