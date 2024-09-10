# python file to create spatial distribution matrices, SDMs
# SDMs defined in excel file, table Tbl_SDM

# SDM is an intermediate matrix used in pathway calculations
# FF1 x SDM
# SDM is used to translate the columns of FF1 into a different spatial resolution,
#  e.g., from receiving grid cells to receiving LMEs.

# calculation depends on whether properties being 'translated' are intensive or extensive

# import geofileops as gfo

from datetime import datetime

import os
import calc.setup.config as cfg
from calc.inout.read_data import getput_idsnames_fromfiles, read_column_from_file, get_data

from calc.functions.gis import read_shapefile_fromfile, project_shapefile
from calc.functions.gis import calc_area, multiintersect

from calc.inout.files import save_df_or_array, read_or_save_pickle
from calc.inout.files import make_filepathext, add_full_directory
from calc.inout.files import get_filepathext_from_files
from calc.functions.checks import checkstr

from calc.functions import stmat

# functions to check that table entries are not empty,
# and that they don't have empty strings or zeros.
# this is useful only because many of the 'x's are nested table/list calls

def create_all_sdmats(sdmat_store_file=''):
    """
    Run through the table of spatial distribution matrices (tbl_stms);
    for each row, create an sdm, if it doesn't already exist
    cfg.sdms is an empty dict created in config.py
    :param sdmat_store_file: name of file (just name, no directory) to store sd_mats
    for later retrieval.  If empty, do not read/save.  If provided, we're reading or saving.
    :return: no return, but does save file if name provided, and updates cfg.sdmats
    """

    print(f'Function <create_all_sdmats> running...')

    fullfilepath = make_filepathext(file=sdmat_store_file,
                                    path=cfg.dir_stms,
                                    ext='.pkl')

    if sdmat_store_file != '':
        # if sdmat_store_file is a file, load sdmats from there:
        if os.path.isfile(fullfilepath):
            print(f'\tFunction <create_all_sdmats> reading '
                  f'd_sdmats from {fullfilepath}...')
            pickleout = read_or_save_pickle(action='read',
                                            file=fullfilepath)
            cfg.dict_STMats = pickleout[0]

            bln_have = True
            bln_save = False
        else:
            bln_have = False
            bln_save = True

    else:
        bln_have = False
        bln_save = False
        # Not saving the full dict_STMats,
        #   but might still save individuals in the create_sdmats function

    if not bln_have:
        # we create them
        print(f'\tFunction <create_all_sdmats> creating sdmats.')
        for idx, row in cfg.dict_excelTables[cfg.tbl_stms].iterrows():
            if row[cfg.tbl_stms_col_create] == True:  # say "==True" to account for possible blanks
                if cfg.bln_debug: print(f'\t\t...working on SDMat row = <{idx}>.  Will either read existing or create.')
                # assign it to dictionary dict_STMats
                cfg.dict_STMats[idx] = create_stmat(tbl_stm_row=row, manual_save_override=False)
            else:
                if cfg.bln_debug: print(f'\t\t...not requested to create SDMat = {idx}')

        # done looping the rows; save what we've created
        if bln_save:
            read_or_save_pickle(action='save', list_save_vars=[cfg.dict_STMats],
                                file=make_filepathext(file=sdmat_store_file, path=cfg.dir_stms, ext='.pkl'))
    else:
        print(f'\t\t...<create_all_sdmats> done')
        print(f'\t\tcfg.d_sdmats keys: {cfg.dict_STMats.keys()}')


def create_sdmat_from_name(nameID, manual_save=False):
    if nameID in cfg.dict_STMats.keys():
        print(f'\t\tsdmat {nameID} already in cfg.sdmats; returning that...')
        return cfg.dict_STMats[nameID]
    else:

        if nameID in cfg.dict_excelTables[cfg.tbl_stms].index:
            datarow = cfg.dict_excelTables[cfg.tbl_stms].loc[nameID]

            sdmat = create_stmat(datarow=datarow, manual_save_override=manual_save)

            if not nameID in cfg.dict_STMats.keys():
                cfg.dict_STMats[nameID] = sdmat
            return sdmat

        else:
            raise KeyError(f'Function <create_sdmat_from_name> called with '
                           f'name = {nameID} that is not in table of sdmats')


def create_stmat(tbl_stm_row, manual_save_override=False):
    """
    For a given row from table tbl_stms, create the sd_mat, depending on type of calc.
    Calculations may be area-based, or include weighting, as described in the STM table.
    See also called function "create_stm_from_df")
    :param tbl_stm_row: a series extracted from table tbl_stms
    :param manual_save_override: if true, prevent from saving; otherwise, use value in table
    :return: a dataframe sdmat, based on parameters in the table
    """

    print(f'\n\t\tFunction <create_stmat> starting...')


    # If the 'Use_Existing_File' column is not empty and not an empty string,
    #  read the file directly and return it.  (And thus end the function)
    if tbl_stm_row[cfg.tbl_stms_col_useExisting] == True:
        # should check that the file exists... but not necessary because <get_data_from_files> will give an error.

        getfile = tbl_stm_row.name  #this gives the index for the row
        if cfg.bln_debug: print(f'\t\tFunction <create_stmat> is returning existing file= {getfile}')

        stm = get_data(datanameID=getfile, return_array=False, fix_ids=False)  #this function should return None if file doesn't exist.
        # could have some better checks that the stm is retrieved properly (e.g., a directory error could make this nothing)

        if not stm is None:
            return stm

    # the file didn't already exist; we need to create it.

    # some text values used only in this function, for renaming columns
    s_orig_ids = 'orig ids'
    s_orig_area = 'area orig'
    s_new_ids = 'new ids'
    s_new_area = 'area new'
    s_weight_ids = 'weight ids'
    s_weight_vals = 'weight values'
    s_weight_area = 'area weight'
    s_isect_area = 'area isect'

    # determine the 'shape': whether the matrix is (orig x new) or vice versa
    if tbl_stm_row[cfg.tbl_stms_col_STMtype] == cfg.str_STMtype_receive:
        stm_shape = stmat.s_stmShape_origXnew
        stm_normalization = stmat.s_stmNormalize_rows
    elif tbl_stm_row[cfg.tbl_stms_col_STMtype] == cfg.str_STMtype_emit:
        stm_shape = stmat.s_stmShape_newXorig
        stm_normalization = stmat.s_stmNormalize_cols
    else:
        raise ValueError('the STM table has the wrong kind of Type', tbl_stm_row[cfg.tbl_stms_col_STMtype])

    orig_geofile = tbl_stm_row[cfg.tbl_stms_col_geoOrig]
    new_geofile = tbl_stm_row[cfg.tbl_stms_col_geoNew]

    # determine if we have weighting (we do if the weight column is not empty)
    bln_have_weight = checkstr(tbl_stm_row[cfg.tbl_stms_col_weightDataUnique])
    if bln_have_weight:
        weight_geofile = cfg.df_data.loc[tbl_stm_row[cfg.tbl_stms_col_weightDataUnique], cfg.tbl_data_col_assocGeo]
        weight_type = cfg.df_data.loc[tbl_stm_row[cfg.tbl_stms_col_weightDataUnique], cfg.tbl_data_col_weightType]
        weight_value_colName = cfg.df_data.loc[tbl_stm_row[cfg.tbl_stms_col_weightDataUnique], cfg.tbl_data_col_headerData]
        weight_ID_colName = cfg.df_files.loc[weight_geofile, cfg.tbl_files_col_headerID]
        # go ahead and get the weight values, which will likely come from a geofile (so we might open twice), but they might not.
        fname, fpath, fext = get_filepathext_from_files(weight_geofile)
        # get the weights -- and don't need to sort them, because we'll merge this onto another frame later.
        ser_weight = read_column_from_file(path=fpath, file=fname, extension=fext,
                                          col=weight_value_colName, colID=weight_ID_colName)
        # we rename the value, but not the ID, since we'll use the IDs to join later
        ser_weight.rename(s_weight_vals, inplace=True) #columns={weight_value_colName: s_weight_vals}
        try:
            ser_weight.index = ser_weight.index.astype('int')
        except ValueError:
            pass
    # determine if we are going to save
    if manual_save_override:
        bln_save = True
    else:
        bln_save = tbl_stm_row[cfg.tbl_stms_col_create]

    # rows and column geofile(s) must be specified, read them
    if cfg.bln_debug2:
        print(f'\t\tFunction <create_stmat> is getting ids, and may need to read a geofile. ')
    master_ids_orig, _ = getput_idsnames_fromfiles(filenameID=orig_geofile, table=cfg.tbl_files)
    master_ids_new, _ = getput_idsnames_fromfiles(filenameID=new_geofile, table=cfg.tbl_files)

    # We're operating on the geofile.
    # We do an intersection to create the df3cols

    list_geofiles = [orig_geofile, new_geofile]
    list_ids_rename = [s_orig_ids, s_new_ids]
    list_areas_rename = [s_orig_area, s_new_area]  # strings used in this function only

    if bln_have_weight:
        list_geofiles.append(weight_geofile)
        list_ids_rename.append(s_weight_ids)
        list_areas_rename.append(s_weight_area)

    list_cols_to_keep = list_ids_rename + list_areas_rename + [s_isect_area]

    if bln_have_weight:
        list_cols_to_keep.append(s_weight_vals)

    temp_isect_name = cfg.str_fileText_Orig + '-' + tbl_stm_row[cfg.tbl_stms_col_geoOrig] + '_' + cfg.str_fileText_New + '-' + tbl_stm_row[cfg.tbl_stms_col_geoNew]
    if bln_have_weight:
        temp_isect_name = temp_isect_name + '_' + cfg.str_fileText_Weight + '-' + weight_geofile

    # check if intersect exists; if so, skip the time-consuming multi-intersect.
    bln_have_isect = os.path.isfile(make_filepathext(file=temp_isect_name, path=cfg.dir_intersects, ext='.pkl'))
    if bln_have_isect:
        shp_intersected = read_or_save_pickle(action='read', file=temp_isect_name, path=cfg.dir_intersects, ext='.pkl')
    else:
        # need to create
        shps = []
        # some of the geographic coordinates can lie exactly upon one another
        #   (e.g., GeosCHEM grid and NEWS2 basins have some line segments shared),
        #   which causes problems for the intersect (creating gaps)
        #   Therefore, apply a small buffer, as per
        #   https://gis.stackexchange.com/questions/277334/shapely-polygon-union-results-in-strange-artifacts-of-tiny-non-overlapping-area/277342
        #   https://gis.stackexchange.com/questions/120286/removing-small-polygon-gaps-in-shapely-polygon/120314#120314
        #   Note that I found it was better to apply the eps buffer only
        #   (rather than remove it with a -eps)
        eps = 0.01

        for i in range(0, len(list_ids_rename)):
            # get file
            if cfg.bln_debug2:
                print(f'\t\tFunction <create_stmat> is reading geofiles for intersection')

            tempshp = read_shapefile_fromfile(geofilenameID=list_geofiles[i], new_id_col=list_ids_rename[i] )

            # convert to the ids to integer, if we can
            try:
                tempshp[list_ids_rename[i]] = tempshp[list_ids_rename[i]].astype(int)
            except ValueError:
                # if we cannot convert to int, we get a value error
                pass

            if cfg.bln_debug2:
                tempshp.to_file(fr'C:\temp\stm\tempshp_noproj_{i}.gpkg', driver='GPKG')
                # tempshp.to_file(fr'C:\temp\stm\shp\tempshp_noproj_{i}.shp',
                #                 driver='ESRI Shapefile')
            # project it
            tempshp = project_shapefile(shp=tempshp, projection_type=cfg.proj_crs_default, projection_string=cfg.proj_s_default)

            if cfg.bln_debug: print(f'\t\tFunction <create_stmat> is buffering')
            # tempshp['geometry'] = tempshp.buffer(eps)  # .buffer(-eps)
            if cfg.bln_debug: print(f'\t\t\t... done.')

            # calc area
            calc_area(shp=tempshp, new_area_name=list_areas_rename[i], conversion_factor=cfg.proj_conv_default)

            if cfg.bln_debug2:
                tempshp.to_file(fr'C:\temp\stm\tempshp_projarea_{i}.gpkg', driver='GPKG')

            # if cfg.bln_debug:
            #     import matplotlib.pyplot as plt
            #     tempshp.plot()
            #     plt.show()

            shps.append(tempshp)


        # what does geofileops do?
        # https://geofileops.readthedocs.io/en/stable/user_guide.html

        # gfo.intersection(input1_path=fr'C:\temp\stm\tempshp_projarea_0.gpkg',
        #               input2_path=fr'C:\temp\stm\tempshp_projarea_1.gpkg',
        #               output_path=fr'C:\temp\stm\tempshp_projarea_0-1.gpkg')
        #
        # gfo.intersection(input1_path=fr'C:\temp\stm\tempshp_projarea_0-1.gpkg',
        #               input2_path=fr'C:\temp\stm\tempshp_projarea_2.gpkg',
        #               output_path=fr'C:\temp\stm\tempshp_projarea_0-1-2.gpkg')


        # run the intersection of the geo... #union to keep all areas
        if cfg.bln_debug:
            print(f'\t---\n\tGoing to run multiintersect with {list_geofiles}')
        # TODO: if shapefiles are the same... avoid intersect?
        shp_intersected = multiintersect(list_shapes=shps, how='union',  # was union  # 5 mins with 'intersection'
                                         new_area_col=s_isect_area,
                                         new_area_conversion=cfg.proj_conv_default)

        # done with multi-intersect...

        # possibly save
        if cfg.bln_save_intersects:
            # print(f'\t\t\t... saving intersect <{temp_isect_name}> as gpkg ...')
            # shp_intersected.to_file(make_filepathext(file= temp_isect_name, path = cfg.dir_intersects, ext= '.gpkg'),
            #                         driver='GPKG')
            print(f'\t\t\t\t... done. Now saving intersect <{temp_isect_name}> pickle...')
            read_or_save_pickle(action='save',
                                file=make_filepathext(file= temp_isect_name, path = cfg.dir_intersects, ext= '.pkl'),
                                list_save_vars=shp_intersected)
            print(f'\t\t\t\t... done.')

    # join the weight values
    if bln_have_weight:
        shp_intersected = shp_intersected.merge(ser_weight, how='left', left_on=s_weight_ids, right_index=True)

    shp_isect_small = shp_intersected.loc[:, list_cols_to_keep]

    #create the stm
    if bln_have_weight:
        stm = stmat.create_stm_from_df(dframe=shp_isect_small,
                                       stm_shape=stm_shape, stm_normalize=stm_normalization,
                                       col_orig_IDs=s_orig_ids, col_new_IDs=s_new_ids,  col_isect_area=s_isect_area,
                                       col_weight_IDs=s_weight_ids,
                                       col_weight_vals=s_weight_vals,
                                       weight_type=weight_type, use_area_weight_if_needed=cfg.bln_area_weight_if_needed)
    else:
        # without weight
        stm = stmat.create_stm_from_df(dframe=shp_isect_small,
                                       stm_shape=stm_shape, stm_normalize=stm_normalization,
                                       col_orig_IDs=s_orig_ids, col_new_IDs=s_new_ids, col_isect_area=s_isect_area)

    # reindex both - even though we reindexed in create_stm_from_df, that function may not have seen *all* the IDs (depends on the intersect)
    # Therefore, we need to reindex to the master lists.
    if stm_shape == stmat.s_stmShape_origXnew:
        stm = stm.reindex(index=master_ids_orig, fill_value=0)
        stm = stm.reindex(columns=master_ids_new, fill_value=0)
    elif stm_shape == stmat.s_stmShape_newXorig:
        stm = stm.reindex(index=master_ids_new, fill_value=0)
        stm = stm.reindex(columns=master_ids_orig, fill_value=0)
    else:
        # can't get here; we alerady trapped for this above
        pass

    # because the geofile indices may not be sorted...
    # we have to sort, whether we create from geofiles or data/files
    # https://stackoverflow.com/questions/17315881/how-can-i-check-if-a-pandas-dataframes-index-is-sorted
    # TODO - could add a check to see if index is sorted
    stm.sort_index(axis='index', inplace=True)
    stm.sort_index(axis='columns', inplace=True)

    if bln_save:
        if checkstr(tbl_stm_row[cfg.tbl_stms_col_fileDirMain]):
            temppath = add_full_directory(tbl_stm_row[cfg.tbl_stms_col_fileDirMain])
            if temppath[-1] != os.sep:
                temppath = temppath + os.sep
        else:
            temppath = cfg.dir_stms

        read_or_save_pickle('save', file=tbl_stm_row[cfg.tbl_stms_col_fileName],
                            path=temppath, ext='.pkl', list_save_vars=stm)
        # save_df_or_array(data=stm, path=temppath,
        #                  filename=tbl_stm_row[cfg.tbl_stms_col_fileName],
        #                  extension=tbl_stm_row[cfg.tbl_stms_col_fileExt])
        if cfg.bln_debug:
            save_df_or_array(data=stm, path = temppath,
                             filename=tbl_stm_row[cfg.tbl_stms_col_fileName],
                             extension='.csv')
    return stm
