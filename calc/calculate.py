import os
from calc.setup import config as cfg
from setup import setup

from calc.setup.create_GEP import adjust_GEP_vector
from calc.inout.read_data import (read_column_from_file, read_matrix_from_file,
                                  getput_idsnames_fromfiles, read_or_save_pickle)
from calc.inout.files import get_filepathext_from_data, write_df_to_excel
from calc.setup.create_stmats import create_stmat

from datetime import datetime
import pandas as pd
import numpy as np
pd.options.mode.copy_on_write = True #https://pandas.pydata.org/pandas-docs/stable/user_guide/copy_on_write.html#migrating-to-copy-on-write

# suffix to attach to output file
cfg.suffix_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

# where are the files: both input and output
# NOTE: also check file locations in config.py
file_settings = r'ozone STM calc setup v00a.xlsx'
file_output = r'ozone STM calc setup v00a'  # will be written as excel file

# Note: which scales to use for impacts (regional, global, or both) are set in excel

# boolean to adjust emission side of FF:
# if True, then use the STM_emit to change rows of FF_raw, to create FF_e
# if False, set FF_e = FF_raw
# bln_adjust_FFe = False
# TODO: include this toggle for both emit and receive, and record this in the output dataframe or filename

# settings to possibly save full matrix files #TODO
bln_save_intermed_FFe = False
bln_save_intermed_FFer = False
bln_save_intermed_FFerxRF = False
bln_save_intermed_FFerxRFxEF = False
bln_save_intermed_FFerxRFxEFxGEPnorm = False

# settings to save output vectors # TODO
bln_save_output_CF = True
bln_save_output_CFxGEPnorm = True

# some settings from the config file.
print(f'\nGlobal parameters: ')
print(f'Area weight if needed: {cfg.bln_area_weight_if_needed}')
print(f'GEP aggregation: {cfg.str_GEPsetup_selected}')
print(f'')

# name of file use to retrieve (or save, if it's missing) the STMs, spatial transform matrices
file_stmats = 'dict_of_stmats'

# ------- Begin calculations -----------

# read the excel file and get all the tables.
setup.main_setup(main_excel_file=file_settings,
                 stmat_file=file_stmats, manual_stmat_skip=False)

if cfg.bln_debug2:
    # look at some of the tables we've created
    # print(cfg.df_stms)
    # print(cfg.dict_excelTables[cfg.tbl_files].head())
    print(cfg.df_calc.head())

# test creation of STMs:
# test1 = create_stmat(datarow =cfg.df_stms.loc['STM newReceive geos to ecoreg', :], manual_save_override=False)
#
# list_stms = ['STM newReceive geos to ecoreg', 'STMw newEmit GeosGrid Countries WtPop', 'STMw newEmit GeosGrid Countries WtAgriNH3',
#              'STMw newEmit GeosGrid Countries WtAgriNOx', 'STMw newEmit GeosGrid Countries WtAgriSO2', 'STMw newEmit GeosGrid Countries WtNoagNH3',
#              'STMw newEmit GeosGrid Countries WtNoagNOx', 'STMw newEmit GeosGrid Countries WtNoagSO2', 'STMw newEmit GeosGrid World WtPop',
#              'STMw newEmit GeosGrid World WtAgriNH3', 'STMw newEmit GeosGrid World WtAgriNOx', 'STMw newEmit GeosGrid World WtAgriSO2',
#              'STMw newEmit GeosGrid World WtNoagNH3', 'STMw newEmit GeosGrid World WtNoagNOx', 'STMw newEmit GeosGrid World WtNoagSO2'
#              ]
#
#
# list_stms = ['STMw newEmit GeosGrid World WtGeneralNH3','STMw newEmit GeosGrid World WtGeneralNOx','STMw newEmit GeosGrid World WtGeneralSO2']
#
# for i in range(0, len(list_stms)):
#     temp = create_stmat(tbl_stm_row=cfg.df_stms.loc[list_stms[i],:])

print(f'\n----------\nAll done with setup.  Starting main calculations')

#
# def prepend_string(value, prefix='id_'):
#     return f"{prefix}{value}"
#
#
# # region ------ Loop the calculation table.  Get the relevant FF, etc., multiply them, and record the results into a overall output dataframe
# rowNum = 0; calcNum = 0
# for calc_idx, calc_row in cfg.df_calc.iterrows():
#     rowNum += 1
#     if calc_row[cfg.tbl_calc_col_doCalc] != True:
#         if cfg.bln_debug:
#             print(f'Not performing main calc table # {rowNum}; the calc column is set to false')
#     else:
#         calcNum += 1
#         print(f'\n----------------------\nPerforming main calc table row number {rowNum}, calc number {calcNum}:')
#         print(f'\t flow = {calc_row[cfg.tbl_calc_col_substance]}; emit = {calc_row[cfg.tbl_calc_col_emitResNew]}; receive = {calc_row[cfg.tbl_calc_col_recvResNew]}')
#
#         # the new emit resolution is important: it sets the number of rows in the CF.  Therefore:
#         # we track each different emission resolution in a separate df, which are held in a dictionary, which we can combine in the end
#         newEmitRes = calc_row[cfg.tbl_calc_col_emitResNew]
#
#         if newEmitRes in cfg.dict_output_by_emitRes.keys():
#             pass
#         else:
#             cfg.dict_output_by_emitRes[newEmitRes] = cfg.df_out_multiindex.copy()  # new blank multiindex dataframe; the emission rows will set the number of rows here.
#
#         # temporary holder, just to make reading easier: dfo = 'data frame output'
#         dfo = cfg.dict_output_by_emitRes[newEmitRes]
#
#         # get IDs for this calculation
#         if cfg.bln_debug: print(f'getting relevant emission IDs and receiving IDs.  '
#                                 f'Will use existing dictionary if we have already read; may need to read some files. ')
#         emission_IDs, emission_Names = getput_idsnames_fromfiles(calc_row[cfg.tbl_calc_col_emitResNew],
#                                                                  table=cfg.tbl_files)
#
#         #TODO: what do we do when there are no names? For now, turn the IDS into names
#         if isinstance(emission_Names, str): #TODO need to think through the logic; currently, returns 'No names for geofile' when there are no names
#             emission_Names = emission_IDs.to_frame()
#             emission_Names[cfg.s_authoritative_names] = emission_Names[cfg.s_authoritative_ids].apply(prepend_string)
#             emission_Names[cfg.s_authoritative_ISO3] = emission_Names[cfg.s_authoritative_names]
#
#         receive_IDs_raw, receive_Names = getput_idsnames_fromfiles(calc_row[cfg.tbl_calc_col_recvResNew],
#                                                                    table=cfg.tbl_files)
#
#         # region -------  Get raw data -------------------
#         # the 'raw' values are read directly from the files, without manipulation
#         # FF, STM come in as matrices and remain matrices
#         # RF, EF come in as vectors, and remain vectors
#         # GEP comes in s vector and is transformed to matrix.  Therefore, this one has "_vector" and "_matrix" in the variable names
#
#         # --- Fate Fators -----
#         # for terrestrial acidifcation:
#         #   all of our FFs are matrices; we could expand this code to get columns and convert (diagonalize) to matrices.
#         #   we also know the FFs do not have any NaNs, so we don't have to check for them.
#         id_FateF = calc_row[cfg.tbl_calc_col_FF]
#         file, path, ext = get_filepathext_from_data(id_FateF)
#         datarow = cfg.df_data.loc[id_FateF, :]
#         FF_raw = read_matrix_from_file(path=path, file=file, extension=ext,
#                                        matrix_has_ids=datarow[cfg.tbl_data_col_matHasIDs],
#                                        geofile=datarow[cfg.tbl_data_col_assocGeo],
#                                        matlabName=datarow[cfg.tbl_data_col_matlabName],
#                                        xlsheet=datarow[cfg.tbl_data_col_excelSheetName],
#                                        return_array=True)
#
#         # ---- Spatial Transform Matrices, STMs ---------
#         # for terrestrial acidification:
#         #   We have created these, and thus we know they do not have any NaNs
#
#         if calc_row[cfg.tbl_calc_col_emitResOrig] == calc_row[cfg.tbl_calc_col_emitResNew]:
#             # we are not adjusting the emission resolution, and so the STM_emit_raw will not be used
#             bln_adjust_FF_emit = False
#             id_STM_emit = 'not used'
#             STM_emit_raw = 0
#
#         else:
#             # get the emit STM
#             #  Note for AH: is it better to do these with function calc.inout.read_data.get_data... which can remove the IDs?
#             # they come in as dataframes, so convert to numpy.  We don't have to track the IDs, since we know they match the authoritative ID lists.
#             bln_adjust_FF_emit = True
#             id_STM_emit = calc_row[cfg.tbl_calc_col_emitSTM]
#             # STM_emit_raw = cfg.dict_STMats[id_STM_emit].to_numpy()
#             STM_emit_raw = cfg.dict_STMats[id_STM_emit]
#
#         # get the receive STM
#         id_STM_recv = calc_row[cfg.tbl_calc_col_recvSTM]
#         # STM_recv_raw = cfg.dict_STMats[id_STM_recv].to_numpy()
#         STM_recv_raw = cfg.dict_STMats[id_STM_recv]
#
#         # ------- Response Factors ---------
#         # as 'raw', they could have NaNs.
#         id_RF = calc_row[cfg.tbl_calc_col_RF]
#         if id_RF in cfg.dict_RFs.keys():
#             RF_raw = cfg.dict_RFs[id_RF]
#         else:
#             file, path, ext = get_filepathext_from_data(id_RF)
#             datarow = cfg.df_data.loc[id_RF, :]
#             RF_raw = read_column_from_file(path=path,file=file, extension=ext,
#                                            col=datarow[cfg.tbl_data_col_headerData],
#                                            sheet=datarow[cfg.tbl_data_col_excelSheetName],
#                                            colID=datarow[cfg.tbl_data_col_headerID],
#                                            fix_IDs=True, return_array=True)
#             cfg.dict_RFs[id_RF] = RF_raw
#
#         # ---- Effect Factors -------
#         # as 'raw', they may have NaNs
#         id_EF = calc_row[cfg.tbl_calc_col_EF]
#         if id_EF in cfg.dict_EFs.keys():
#             EF_raw = cfg.dict_EFs[id_EF]
#         else:
#             file, path, ext = get_filepathext_from_data(id_EF)
#             datarow = cfg.df_data.loc[id_EF, :]
#             EF_raw = read_column_from_file(path=path,file=file, extension=ext,
#                                            col=datarow[cfg.tbl_data_col_headerData],
#                                            sheet=datarow[cfg.tbl_data_col_excelSheetName],
#                                            colID=datarow[cfg.tbl_data_col_headerID],
#                                            fix_IDs=True, return_array=True)
#             cfg.dict_EFs[id_EF] = EF_raw
#
#         # ------ GEP ------------
#         # as 'raw', they may have NaNs
#         id_GEP = calc_row[cfg.tbl_calc_col_GEP]
#         if id_GEP in cfg.dict_GEPs.keys():
#             GEP_raw_vector = cfg.dict_GEPs[id_GEP]
#         else:
#             file, path, ext = get_filepathext_from_data(id_GEP)
#             datarow = cfg.df_data.loc[id_GEP, :]
#             GEP_raw_vector = read_column_from_file(path=path, file=file, extension=ext,
#                                                    col=datarow[cfg.tbl_data_col_headerData],
#                                                    sheet=datarow[cfg.tbl_data_col_excelSheetName],
#                                                    colID=datarow[cfg.tbl_data_col_headerID],
#                                                    fix_IDs=True,
#                                                    return_array=True)
#             cfg.dict_GEPs[id_EF] = GEP_raw_vector
#
#
#         # we have have now gotten raw FF, STM_emit, STM_recv, RF, EF, and GEP
#         # endregion (Get raw data)
#
#         # region ------- Change FF resolution --------
#         # adjust the FF emission and receiving resolution by multiplying by the STMs
#         #
#         if bln_adjust_FF_emit:
#             FF_e_raw = np.matmul(STM_emit_raw, FF_raw)
#         else:
#             FF_e_raw = FF_raw
#
#         # adjust the receiving resolution
#         FF_e_r_raw = np.matmul(FF_e_raw, STM_recv_raw)
#         # endregion
#
#         # region ------ Drop receiving locations with NaNs in RF, EF, and GEP ----------
#         # Check for NaNs (in receiving vectors) and remove these IDs to do the matrix operations.
#         #   We'll later put the IDs back to the authoritative list, which will leave NaN as NaN, rather than 0.
#         #   (i.e., we'll be clear that the value is unknown, not 0).
#
#         combo_df = pd.concat([receive_IDs_raw.to_frame(name='IDs'),
#                               pd.DataFrame(RF_raw, columns=['RF']),
#                               pd.DataFrame(EF_raw, columns=['EF']),
#                               pd.DataFrame(GEP_raw_vector, columns=['GEP']),
#                               ], axis=1)
#
#         if cfg.bln_debug:
#             print(f'Number of NaNs in the raw RF, EF, GEP:\n{combo_df.isna().sum()}')
#
#         # Create boolean mask (as a series) for rows with NaN values
#         recv_nan_mask = combo_df.isnull().any(axis=1)
#
#         FF_e_r = FF_e_r_raw[:, ~recv_nan_mask.values]
#
#         # for the vectors, drop the masked entries (which have an NaN in one of the 3 vectors) and make sure they are 1 x n vectors.
#         RF = RF_raw[~recv_nan_mask.values]
#         if RF.ndim==1:  #this is a numpy series, i.e., 1-D, but without a number of columns.
#             RF = RF.reshape(1,-1)  # convert to a 1 x n array
#         elif RF.ndim==2:  # this has two dimensions; check if the vector is a column vector (n x 1); if so, transpose to row vector (1 x n)
#             if RF.shape[1]==1: RF = RF.T
#
#         EF = EF_raw[~recv_nan_mask.values]
#         if EF.ndim==1:
#             EF = EF.reshape(1,-1)
#         elif EF.ndim==2:
#             if EF.shape[1]==1: EF = EF.T
#
#         GEP_noNaN_vector = GEP_raw_vector[~recv_nan_mask.values]
#         # Normalize GEP and convert to matrix
#         GEP_normalized_matrix = adjust_GEP_vector(series=pd.DataFrame(GEP_noNaN_vector, columns=['GEP']),
#                                                   adjust_method=calc_row[cfg.tbl_calc_col_GEPnorm_method],
#                                                   num_emit_rows=len(emission_IDs))
#
#         # endregion (Drop receiving locations with NaNs in RF, EF, and GEP)
#
#
#         # region ---- FF to CF multiplication --------
#         # perform the matrix multiplication.
#         # Results are now as (emit x receive) matrices, so total CFs are based on summing along the rows to create (n x 1) vectors
#
#         # RF and EF are the Hadamard product: element-by-element.
#         # numpy "*" broadcasts the vector across the matrix.  (Note that the vectors have to be (1 x n) row vectors to be brodcast correctly.
#         FFxRF = FF_e_r * RF
#         FFxRFxEF = FFxRF * EF
#         FFxRFxEFxGEP = FFxRFxEF * GEP_normalized_matrix
#
#         # totals
#         # averaging across ecoregions
#         CF_regional = np.average(FFxRFxEF, axis=1)
#         CF_gep = np.average(FFxRFxEFxGEP, axis=1)
#         # additional step if emission resolution is "World" (calculation of the global average value)
#         if calc_row[cfg.tbl_calc_col_emitResNew] == 'World':
#             countries_IDs, countries_Names = getput_idsnames_fromfiles('Countries',
#                                                                      table=cfg.tbl_files)
#             CF_regional = CF_regional / len(countries_Names)
#             CF_gep = CF_gep / (len(countries_Names)**2)
#
#         # save the CFs into df
#         df_CF = emission_Names  # TODO: IF these are not named... then what?  (e.g., if 'native' CFs with GeosGrid)
#         df_CF['CF ' + cfg.str_impactScale_regional] = CF_regional
#         df_CF['CF ' + cfg.str_impactScale_global] = CF_gep
#
#         # TODO: add saving some of the intermediate matrices.
#         # endregion (FF to CF multiplication)
#
#
#         # region ------ Populate the final output matrix
#
#         # saving the CF into the output matrix needs to be done for both regional and global CFs, so loop over these.
#         for impactScale_idx, impactScale_row in cfg.df_impactScale.iterrows():
#
#             if impactScale_idx == cfg.str_impactScale_regional:
#                 CF = CF_regional
#             elif impactScale_idx == cfg.str_impactScale_global:
#                 CF = CF_gep
#             else:
#                 # Shouldn't be able to reach this point... but somehow, the impact scale did not match the possible options...
#                 CF = None
#
#             if impactScale_row[cfg.tbl_impactScale_col_Include]==True: # thet we are recording this calc scale
#
#                 # Create a multiindex for accessing data in the dfo.
#                 # NOTE: must specify the correct order of multiindex here, and the order is set in config:
#                 #   list_multiindex_col_levels = [multiindex_col_emitResNew, multiindex_col_subst, multiindex_col_impactScale, multiindex_col_sectorWeight]
#                 temp_multiindex = (calc_row[cfg.tbl_calc_col_emitResOrig],  # original emission resolution
#                                    calc_row[cfg.tbl_calc_col_emitResNew],  # new emission (aggregation) Resolution
#                                    calc_row[cfg.tbl_calc_col_substance],  # Substance
#                                    impactScale_idx,  # Impact Scale:
#                                    calc_row[cfg.tbl_calc_col_emitWeightSector])  # sector weighting
#
#                 # record data into the dfo
#                 dfo.loc[:, (*temp_multiindex, cfg.multiindex_col_dataItem_CF)] = CF
#                 # Need to track names as well, because the emission IDs, names, etc. can change with aggregation resolution (globe or country).
#                 #   Thus, emission ID could be ambiguous at different resolutions.  We assume the names are not...
#                 # the CF was added as an array without an index, creating a 0-n index in the receiving dataframe.  I.e., not based on emission IDs.
#                 # to prevent the emission names from matching on their index (which isn't necessarily 0-n), use .values
#                 dfo.loc[:, (*temp_multiindex, cfg.multiindex_col_dataItem_Name)] = emission_Names[cfg.s_authoritative_names].values
#                 dfo.loc[:, (*temp_multiindex, cfg.multiindex_col_dataItem_Code)] = emission_Names[cfg.s_authoritative_ISO3].values
#
#                 if cfg.bln_debug2:
#                     print(dfo.head(5))
#
#         # done looping the impact Scales
#         # now record the dfo back into the dictionary of output, indexed by emitting resolution
#         cfg.dict_output_by_emitRes[newEmitRes] = dfo  #write back into the dictionary, to be accessed again on the next loop of calc_row
#
#         # endregion
#
# # endregion (Loop calculation table)
#
# count_emitResKeys = 0
# for emitResKeys in cfg.dict_output_by_emitRes.keys():
#     count_emitResKeys += 1
#     # get the relevant df and hold in a short variable name.
#     dfo = cfg.dict_output_by_emitRes[emitResKeys]
#
#     # region --- calculate general sector weighting
#     # if requested go back through and calculate the general sector:
#     # can only do this if both ag and non-ag have been calculated:
#     list_sectors_used = dfo.columns.unique(level=cfg.multiindex_col_sectorWeight).to_list()
#     if (cfg.str_sector_agric in list_sectors_used) and (cfg.str_sector_notAg in list_sectors_used):
#
#         for emit_res_orig in dfo.columns.unique(cfg.multiindex_col_emitResOrig).to_list():
#
#             for emit_res_new in dfo.columns.unique(cfg.multiindex_col_emitResNew).to_list():
#
#                 for subst in dfo.columns.unique(cfg.multiindex_col_subst).to_list():
#
#                     for impScale in dfo.columns.unique(cfg.multiindex_col_impactScale).to_list():
#
#                         # check if both exist for this combination:
#                         if (dfo.columns.isin([(emit_res_new, subst, impScale, cfg.str_sector_agric)]).any()) and (
#                         dfo.columns.isin([(emit_res_new, subst, impScale, cfg.str_sector_notAg)]).any()):
#
#                             if cfg.bln_debug:
#                                 print(f'\t\tCalculating general weighting; we can do so for {emit_res_new, subst, impScale}')
#
#                             temp_agric = dfo.loc[:, (emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_agric, cfg.multiindex_col_dataItem_CF)]
#                             temp_nonag = dfo.loc[:, (emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_notAg, cfg.multiindex_col_dataItem_CF)]
#
#                             temp_genrl = (temp_agric + temp_nonag) / 2  #average the ag and non-ag
#
#                             dfo.loc[:,(emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_genrl, cfg.multiindex_col_dataItem_CF)] = temp_genrl
#                             # also record the names and IDs.  Can take from either agric or nonagric.
#                             dfo.loc[:, (emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_genrl, cfg.multiindex_col_dataItem_Name)] = \
#                                 dfo.loc[:, (emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_agric, cfg.multiindex_col_dataItem_Name)]
#                             dfo.loc[:, (emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_genrl, cfg.multiindex_col_dataItem_Code)] = \
#                                 dfo.loc[:, (emit_res_orig, emit_res_new, subst, impScale, cfg.str_sector_agric, cfg.multiindex_col_dataItem_Code)]
#                             if cfg.bln_debug2:
#                                 print(f'\t... new general: ')
#                                 print(dfo.columns)
#
#         # so everything is grouped by substance, etc., sort the columns and the data by the indices
#         dfo = dfo.sort_index(axis=1, level=[0, 1, 2, 3, 4, 5], ascending=[True, True, True, True, True, True])
#
#         if cfg.bln_debug:
#             print(f'\tAfter adding general columns and sorting, here is the column list:')
#             print(dfo.columns)
#     else:
#         # agric and non-agric not present, don't both looking to do the averaging
#         pass
#
#     # now DFO has a 6-level multiindex columns. The last level, 'Data Item' has CFs, location names, and location IDs.
#
#     # endregion
#
#     # region ----- Rearrange and add output
#
#     # stack brings several index levels at the column to a row index.
#     # NOTE: check we get all levels of multiindex here.  (If not, they appear as levels in the column indices)
#     df_stack = dfo.stack(level=[cfg.multiindex_col_emitResOrig, cfg.multiindex_col_emitResNew, cfg.multiindex_col_subst, cfg.multiindex_col_impactScale, cfg.multiindex_col_sectorWeight], future_stack=True)
#
#     # because the row index already existed, we now have an unneeded multi-index on the rows.
#     # the row index was just a 0-n, so we can drop it -- it wasn't tied to emission locations.
#     df_stack.index = df_stack.index.droplevel(level = 0)
#
#     # easier to move everything to columnns, where it can be used in the merge.  This will create another dummy index, 0-n.
#     df_build_unsort = df_stack.reset_index(level=list(range(0, df_stack.index.nlevels)))
#
#     # to match the typical layout, sort so that we see the lists of emitting areas:
#     df_build = df_build_unsort.sort_values([cfg.multiindex_col_subst, cfg.multiindex_col_impactScale, cfg.multiindex_col_sectorWeight], ascending = [True, True, True])
#
#     # so index is nice... starts from 0
#     df_build = df_build.reset_index(drop=True)
#
#     # merge data in; write values
#     # merge flow name and CAS: will merge on the substance
#     df_build = pd.merge(left=df_build, right=cfg.df_flows.loc[:, [cfg.str_outputCol_FlowName, cfg.str_outputCol_FLowCAS]],
#                         how='left', left_on = cfg.multiindex_col_subst, right_index=True)
#
#     # merge units
#     df_build = pd.merge(left=df_build, right=cfg.df_impactScale.loc[:, cfg.str_outputCol_unit], left_on=[cfg.multiindex_col_impactScale], right_index=True)
#
#     df_build[cfg.str_outputCol_CF_uncert_low] = ''
#     df_build[cfg.str_outputCol_CF_uncert_high] = ''
#
#     df_build[cfg.str_outputCol_flowclass0] = ''
#     df_build[cfg.str_outputCol_flowclass1] = ''
#     df_build[cfg.str_outputCol_flowclass2] = ''
#
#     df_build[cfg.str_outputCol_species] = 'Aggregated'
#     df_build[cfg.str_outputCol_lciaMethodRealm] = 'Terrestrial'
#     df_build[cfg.str_outputCol_lciaMethodMathAppr] = 'Marginal'
#
#     df_build[cfg.str_outputCol_scenario] = ''
#     df_build[cfg.str_outputCol_CF_deriv] = ''
#
#     df_build[cfg.str_outputCol_matchCF] = ''
#     df_build[cfg.str_outputCol_matchComp] = ''
#
#     # reorder the columns to match desired output
#     new_col_order1 = cfg.list_out_cols
#     new_col_order2 = [c for c in list(df_build.columns) if c not in cfg.list_out_cols]
#
#     df_emitRes_final = df_build.reindex(columns=new_col_order1 + new_col_order2)
#
#     # endregion (rearrange output)
#
#     if count_emitResKeys == 1:
#         cfg.df_out_final = df_emitRes_final
#     else:
#         # stack the new out_final to the existing
#         cfg.df_out_final = pd.concat([cfg.df_out_final, df_emitRes_final], axis=0)
#
# # region --------- Final output ---------
# str_output_file = file_output + '_' + cfg.suffix_time
# write_df_to_excel(df=cfg.df_out_final, file = str_output_file, path = cfg.dir_output, ext='.xlsx', merge_cells=False)
#
# # info about run
# # TODO: add to this.  Perhaps add a sheet with the used calculation table, etc.
# description_data = [
#     {"Item": "Script", "Description": os.path.basename(__file__)},
#     {"Item": "Run date", "Description": cfg.suffix_time}
# ]
#
# # Create the DataFrame
# df_info = pd.DataFrame(description_data)
#
# write_df_to_excel(df=df_info, file = str_output_file, path = cfg.dir_output, ext = '.xlsx', sheet_name='Calc Info')
#
# df_calc = cfg.df_calc[cfg.df_calc[cfg.tbl_calc_col_doCalc] == True]
# write_df_to_excel(df=df_calc, file = str_output_file, path = cfg.dir_output, ext = '.xlsx', sheet_name='Calc Table')


# endregion

# # create output dictionary
# output_dict = OrderedDict([
#     ('Calculation information --------', ''),
#     ('Date/time', cfg.suffix_time),
#     ('Settings file', file_settings),
#     ('Aggregation --------', ''),
#     ('Is general sector based on average of calculated agric/nonagric?', aggreg_general_from_others),
#     ('Use area-weighting with NaN likelihood?', area_weight_with_nan_likelihood),
#     ('Use area-weighting with 0 or Nan likelihood?', area_weight_with_0orNan_likelihood),
#     ('Aggregate from target to lower resolution? (e.g., from countries to continents)', aggregate_up),
#     ('Freshwater Eutrophication --------', ''),
#     ('Add a soil factor based on freshwater?', add_fw_eutroph_soil),
#     ('If so, factor=', fw_eutroph_soil_factor),
#     ('Absolute/Relative factors --------', ''),
#     ('Calculated relative factors as well?', normalize_results)
#     ]
#     )

