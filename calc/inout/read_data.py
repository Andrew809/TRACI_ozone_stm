import os  # check for file existence, use os.sep
from collections import OrderedDict

import geopandas as gpd

# need pyogrio
# conda install -c conda-forge pyogrio

import numpy as np

import pandas as pd
import scipy
from scipy import io

from calc.setup import config as cfg
from calc.inout.files import (get_filetablecol, get_filepathext_from_files, get_filepathext_from_data,
                              make_filepathext, read_or_save_pickle)
from calc.inout.files import add_full_directory, combine_dir_subdir

from calc.functions.checks import checkstr

def get_fileUniqueID_from_pathFileExt(path, file, extension):
    thisFilePathExt = make_filepathext(file, path, extension)
    # thisFileUniqueID = cfg.df_files.loc[cfg.df_files[cfg.tbl_files_colAdd_pathFileExt] == thisFilePathExt].index.tolist()[0]
    mask = cfg.df_files[cfg.tbl_files_colAdd_pathFileExt] == thisFilePathExt
    fileUniqueID = cfg.df_files.index[mask].tolist()[0]
    return fileUniqueID

def try_convert_ids_to_int(theids):
    """
    convert to the ids to integer, if we can.
    Necessary because some IDs are stored as strings...

    :param theids:
    :return:
    """
    try:
        theids = theids.astype('int')
    except ValueError:
        # if we cannot convert to int, we get a value error
        pass
    return theids

def get_data(datanameID, return_array=False, fix_ids=True, fix_nans=False):
    """
    # given a datanameID , pull out the data from tbl_stms or tbl_data
    # returning as a numpy array if so specified
    :param datanameID: the row index in the table of data
    :param return_array: boolean; if true, return as numpy array
    :return: a series (with ids) if a column, or numpy array (no ids) of
    the data associated with datanameID
    """
    print(f'\tFunction <get_data> called datanameID = <{datanameID}> and with return_array={return_array}')

    # get_filetablecol tells us the corresponding source file nameID and the
    #  appropriate table (either sheetfiles or geofiles)
    # e.g., for a Tbl_Data entry of FFs, we need to know which excel file to look in
    filenameID, table, col = get_filetablecol(datanameID=datanameID)

    if filenameID is None:
        raise KeyError(f'In function <get_data>, datanameID = {datanameID} not found in tables that <get_filetablecol> searches')

    else:
        # we know it's in either stms or data
        # in tbl stms ----------
        if table == cfg.tbl_stms:
            filename = cfg.df_stms.loc[filenameID, cfg.tbl_stms_col_fileName]
            if checkstr(cfg.df_stms.loc[filenameID, cfg.tbl_stms_col_fileDirMain]):
                path = combine_dir_subdir(cfg.df_stms.loc[filenameID, cfg.tbl_stms_col_fileDirMain],
                                          cfg.df_stms.loc[filenameID, cfg.tbl_stms_col_fileDirSub])
            else:
                path = cfg.dir_stms
                if path[-1] != os.sep:  path += os.sep  #add a trailing path separator
            ext = cfg.df_stms.loc[filenameID, cfg.tbl_stms_col_fileExt]

            # we know we created this stm, so we can read it without much checking
            if os.path.isfile(path+ os.sep + filename + ext):
                the_stm = read_matrix_from_file(path=path, file=filename,
                                                extension=ext,
                                                matrix_has_ids=True,
                                                return_array=True)
                return the_stm
            else:
                return None

        # in tbl data -----------------
        if table == cfg.tbl_files:
            filename, path, ext = get_filepathext_from_data(datanameID=datanameID)

            if not os.path.isfile(path + filename + ext):
                raise FileExistsError(f'Function <get_data> did not find file {path + filename + ext}, '
                                      f'\n\tassociated with {datanameID} in datatable and {filenameID} in table {table}')


            row_tbl_files = cfg.df_files.loc[filenameID, :]
            row_tbl_data = cfg.df_data.loc[datanameID, :]

            associatedGeoFile = row_tbl_data[cfg.tbl_data_col_assocGeo]

            # check if is matrix
            # to handle cell being blank or not being set to true...

            if row_tbl_data[cfg.tbl_data_col_isMatrix]:

                matlabName = row_tbl_data[cfg.tbl_data_col_matlabName]
                hasids = row_tbl_data[cfg.tbl_data_col_matHasIDs]

                return read_matrix_from_file(path=path, file=filename,
                                             extension=ext,
                                             matrix_has_ids=hasids,
                                             geofile=associatedGeoFile,
                                             matlabName=matlabName,
                                             return_array=return_array)

            # there is no matrix, so this is a standard file retrieval
            thesheet = row_tbl_data[cfg.tbl_data_col_excelSheetName]
            theColID = row_tbl_data[cfg.tbl_data_col_headerID]

            # but we may need to fix the ids, because many datafiles may lack data for some geo ids.
            # So get the preliminary data

            tempdata = read_column_from_file(path=path,
                                             file=filename,
                                             extension=ext,
                                             sheet=thesheet,
                                             col=col,
                                             colID=theColID)

            # sometimes there are n/a, or NA(), etc. col_values in excel files, so we replace with zero
            if fix_nans:
                tempdata.fillna(0, inplace=True)

            if fix_ids:   #then send it through the id fixer
                tempdata = fix_data_ids(df=tempdata,
                                        ref_geofile=associatedGeoFile,
                                        sort=True)

            if return_array:
                return tempdata.to_numpy()
            else:
                return tempdata


# def get_data_fromfiles(filenameID, return_array=False):
#     """
#     reads an entire excel sheet or matrix.
#     HasIDs is outside of the 'is matrix' check,
#     because we may want to assume IDs are on an excel sheet.
#     :param filenameID: name of file (row) from ONLY the regular files table
#     :param return_array: boolean; return numpy array if true
#     :return: data as column or matrix, as called.
#     """
#     # we know table is either cfg.t_reg_files
#     # this is for reading the entire sheet or entire matrix.
#
#     if filenameID in cfg.df_files.index:
#         # get file inpfo so we can read
#         filename, path, ext = get_filepathext_from_files(filenameID=filenameID)
#
#         if not os.path.isfile(path + filename + ext):
#             raise FileExistsError(f'Function <get_data> '
#                                   f'did not find file {path + filename + ext}, '
#                                   f'\n\tassociated with {filenameID} '
#                                   f'in table {cfg.tbl_files}')
#
#         tblrow = cfg.df_files.loc[filenameID, :]
#         hasids = tblrow[cfg.tbl_data]
#
#         # to handle cell being blank or not being set to true...
#         if not tblrow[cfg.s_regfile_ismatrix] is None and \
#                 tblrow[cfg.s_regfile_ismatrix]:
#
#             # this is a matrix
#             # check if there associated geo file
#             associatedGeoFile = tblrow[cfg.s_regfile_assocgeo]
#             matlabName = tblrow[cfg.s_regfile_matname]
#
#             return read_matrix_from_file(path=path, file=filename,
#                                          extension=ext,
#                                          matrix_has_ids=hasids,
#                                          geofile=associatedGeoFile,
#                                          matlabName=matlabName,
#                                          return_array=return_array)
#
#         else:
#             thesheet = tblrow[cfg.s_regfile_xlsheet]
#
#             return read_matrix_from_file(path=path, file=filename,
#                                          extension=ext,
#                                          matrix_has_ids=hasids,
#                                          xlsheet=thesheet,
#                                          return_array=return_array)
#
#     else:
#         raise KeyError(f'In function <get_data_fromfiles>, '
#                        f'\n\tfile name = {filenameID} not found in '
#                        f'index of table {cfg.t_reg_files}')


def fix_data_ids(df, ref_geofile, sort=True):
    """
    based on ids in reference geofile, reindex and sort the dataframe
    :param df: dataframe to be reindexed.  (or series)
    :param ref_geofile: name of geofile (in the table of geofiles) with the 'authoritative' list of ids
    :param sort: boolean: do we sort the indices
    :return: reindexed, and possibly sorted, dataframe
    """

    ids, _ = getput_idsnames_fromfiles(filenameID=ref_geofile,
                                       table=cfg.tbl_files,
                                       return_dict=False)

    # determine whether we also need to reindex columns, which should only be the case for a square matrix
    if isinstance(df, pd.Series):
        bln_fix_cols = False
    elif isinstance(df, pd.DataFrame):
        # assume that we can reindex a square matrix...
        # this could be problematic in cases with few IDs and lots of other columns,
        # if they have to be the same number
        bln_fix_cols = df.shape[0] == df.shape[1]
    else:
        bln_fix_cols = False
        raise TypeError(f'Function <fix_data_ids> called with a '
                        f'data type ({type(df)})'
                        f'that is not pd.Series or pd.DataFrame.')

    df2 = df.reindex(index=ids, fill_value=0)

    if bln_fix_cols:
        df2 = df2.reindex(columns=ids, fill_value=0)

    if sort:
        df2.sort_index(inplace=True)
        if bln_fix_cols:
            df2.sort_index(inplace=True)

    return df2


def read_column_from_file(path, file, extension, col, sheet='',
                          colID='', fix_IDs=True, return_array=False, realIDcolumn = ''):
    """
    Read columns from xlsx, csv, mat, npy, or shp files
    We assume excel file structured with column headings.
    If data are a column, try to return a pd Series with index being the IDs

    For non-matrix files, this returns a series

    :param path:
    :param file:
    :param extension:
    :param col: column from which to get excel data; if column is empty,
    we should be reading a matrix from matlab, numpy
    :param sheet: excel sheet; if empty, we have a shapefile
    :param colID: for excel or csv, associated column with IDs corresponding to the data
    :param return_array: boolean; if true, return numpy array
    :return: column data: series with index being the IDs;
    matrix data: dataframe with index/columns being IDs
    """

    colExtract = ''

    if path[-1] != os.sep: path += os.sep  # add a trailing separator, if needed

    # check if we already have read the file
    thisFileUniqueID = get_fileUniqueID_from_pathFileExt(path, file, extension)
    bln_have_read = thisFileUniqueID in cfg.dict_filesHaveBeenRead.keys()

    if extension == '.shp' or extension == '.gpkg' :
        # the ID is returned as part of the geodataframe
        if bln_have_read:
            shp = cfg.dict_filesHaveBeenRead[thisFileUniqueID]
        else:
            if os.path.isfile(make_filepathext(thisFileUniqueID, cfg.dir_input_pkls, '.pkl')):
                if cfg.bln_debug: print(f'\t\tReading <{thisFileUniqueID}> from *pickle* in function <read_shapefile_fromfile>')
                shp = read_or_save_pickle('read', cfg.dir_input_pkls + os.sep + thisFileUniqueID + '.pkl')
                if cfg.bln_debug: print(f'\t\t\t...done reading.')
            else:
                print(f'\t\t<read_column_from_file> is reading shapefile {file + extension} from file...')
                shp = gpd.read_file(path + file + extension, engine='pyogrio')  #TODO: speed up: https://gis.stackexchange.com/questions/469503/speed-up-reading-gpkg-as-geopandas-dataframe
                print('\t\t\t...done reading')

                #kludgily deal with data that doesn't have an ID stored as a column...

                if realIDcolumn not in shp.columns:
                    # happens with a file created from a raster; the geopackage ID may not be in columns.
                    # the first time we read, we want to force this column to be included in the columns, not just as an index

                    shp = shp.reset_index()
                    shp = shp.rename(columns={'index': realIDcolumn})


                if colID != '':
                    if colID not in shp.columns:
                        shp = shp.reset_index()
                        shp = shp.rename(columns={'index': colID})


                # we may save them, either in the run-time dictionary of what we've read, and/or as a pickle of the shapefile
                if cfg.bln_keep_filesHaveBeenRead:
                    cfg.dict_filesHaveBeenRead[thisFileUniqueID] = shp.drop(columns='geometry')
                if cfg.bln_save_geopickles:
                    os.makedirs(cfg.dir_input_pkls, exist_ok=True)
                    read_or_save_pickle('save', cfg.dir_input_pkls + os.sep + thisFileUniqueID + '.pkl',
                                        list_save_vars=shp)
                    # and geopackage for later use
                    shp.to_file(cfg.dir_input_pkls + os.sep + thisFileUniqueID + '.gpkg', driver='GPKG')

        if col in shp.columns:
            # we're extracting from a geodataframe, so we get the series
            if colID == '':
                colExtract = shp.loc[:, col]
            else:

                dfExtract = shp.loc[:, [colID, col]]
                try:
                     dfExtract[colID] = dfExtract[colID].astype('int')  # int64 allows for NaNs; we should not have these in the indices!
                except ValueError:
                    # if we cannot convert to int, we get a value error
                    pass
                dfExtract.set_index(keys=colID, drop=True, inplace=True)
                colExtract = dfExtract.squeeze() # turns the single column dataframe into a serios

        else:
            raise KeyError(f'function read_column_from_file '
                           f'could not find column {col} in '
                           f'shapefile {file} with column list: {shp.columns}')

    elif extension == '.xlsx' or extension == '.csv':
        # use pandas to read
        # excel *can* be problematic.  We might have duplicate IDs :(

        if bln_have_read:
            df = cfg.dict_filesHaveBeenRead[thisFileUniqueID]
        else:

            print(f'\t\treading excel or csv for col= {col}, '
                  f'colID= {colID} (okay if colID blank)')
            if extension == '.xlsx':
                df = pd.read_excel(io=path + file + extension, sheet_name=sheet)
            else: # extension == '.csv':
                df = pd.read_csv(filepath_or_buffer=path + file + extension,
                                  header=0)
            if cfg.bln_keep_filesHaveBeenRead:
                cfg.dict_filesHaveBeenRead[thisFileUniqueID] = df


        if col in df.columns:
            # note that vals is extracted as a series, so it has a (0..n) index already.
            vals = df.loc[:, col]

            # even though colID default is '', if we read an empty ID from
            # an excel sheet, then colID is None, so we check using checkstr
            if checkstr(colID):
                idx = df.loc[:, colID]
                colExtract = pd.Series(data=vals.to_numpy(), index=idx)
                # since we have indices, we sort
                colExtract.sort_index(axis=0,inplace=True)

                if colExtract.index.duplicated().any():
                    print('\tWARNING: Duplicate indices found in excel file {path+file+extension}, '
                          'looking at col ID = {colID}; col data = {col}.')

                    # We try to remove duplicates -- if the data column values are the same, too.
                    # Method to remove duplicates where the duplicate indices have the same values across the rows
                    # Reset index, drop duplicates considering the entire DataFrame, and set the index back
                    df_reset = colExtract.reset_index()
                    df_dedup = df_reset.drop_duplicates().set_index(colID)

                    if df_dedup.index.duplicated().any():
                        raise KeyError ('Duplicate indices but NON-duplicate data found in excel file {path+file+extension}, '
                                        'looking at col ID = {colID}; col data = {col}.')
                    else:
                        print(f'\tWarning averted: duplicates in index were also present in the data, so these rows were removed')
                        colExtract = df_dedup.squeeze()  # after moving the index into the series and back out, we have a dataframe .  Squeeze converts back to series.

            else:
                # we cannot sort by IDs
                colExtract = pd.Series(data=vals)

        else:
            raise KeyError(f'Function <read_column_from_file> '
                           f'could not find column {col} '
                           f'in file = {file+extension} with column list: {df.columns}')

    # have not yet added getting IDs from mat or npy
    # elif extension == '.mat':
    #
    #     loaded = io.loadmat(path + file + extension)
    #     # print(loaded.keys())
    #     mat = loaded[matlabName]  # need to get the right variable out of the mat file
    #     return 'Not yet getting IDs'
    #
    # elif extension == '.npy':
    #     # we do not grab any IDs in this function
    #     npy = np.load(path + file + extension)
    #     return 'Not yet getting IDs'

    else:
        raise KeyError(f'function <read_column_from_file> was '
                       f'given an unknown extension: {extension}')

    # now we have defined the col extract,
    # # and we can return, based on dataframe or array

    #TODO - I have the IDs here... should just fix w/out calling fix_IDs
    if fix_IDs:
        assoc_geofile = cfg.df_files.loc[thisFileUniqueID, cfg.tbl_files_col_assocGeo]
        colExtract = fix_data_ids(df=colExtract, ref_geofile=assoc_geofile, sort=True)

    if not isinstance(colExtract, str):
        if return_array:
            return colExtract.to_numpy()
        else:
            return colExtract


def read_matrix_from_file(path, file, extension, matrix_has_ids=False,
                          geofile='', matlabName='', xlsheet='',
                          return_array=False):

    # return a matrix from the table tregfiles (table regular files)
    # if return_array = True, we return just the numbers
    # if not, we return a dataframe.
    # if matrix_has_ids is true, then the file itself has ids, and we set these as index and column
    # if matrix_has_ids is false, then we grab ids from the geofile
    # if return_array is false (we need ids), but matrix_has_ids=False and geofile='',
    #    then we return an error

    print(f'\tFunction <read_matrix_from_file> called with return_array={return_array}')

    if not return_array and not (matrix_has_ids or geofile != ''):
        raise TypeError(f'Function <read_matrix_from_file> called to return a '
                        f'dataframe, but no ids or geofile supplied.'
                        f'\tfile = "{file+extension}", in directory = {path}')

    mat = ''
    bln_got_df = False

    if path[-1] != os.sep: path += os.sep  #add a trailing separator

    print(f'\t\tLoading "{file + extension}" in function read_matrix_from file...')

    if extension == '.mat':
        loaded = io.loadmat(path + file + extension)
        # print(loaded.keys())
        mat = loaded[matlabName]  # need to get the right variable out of the mat file

    elif extension == '.npy':
        mat = np.load(path + file + extension)

    elif extension == '.csv':
        # mat = pd.read_csv(filepath_or_buffer=path+file+extension,
        #                   sep=',', index_col=None, header=None)
        # mat = mat.col_values
        mat = np.genfromtxt(path + file + extension, delimiter=',')

    elif extension == '.xlsx':
        df = pd.read_excel(path + file + extension, sheet_name=xlsheet)
        bln_got_df = True
        mat = df.to_numpy()

    elif extension == '.pkl':
        mat = read_or_save_pickle(action='read', file=file, path=path, ext=extension)
        if isinstance(mat, pd.DataFrame):
            bln_got_df = True
    else:
        raise TypeError(f'function <read_matrix_from_file> called with unknown'
                        f'extension (not mat, npy, csv, pkl, or xlsx)'
                        f'\nfile = {file + extension} in directory = {path}')

    if isinstance(mat, str) and mat == '':
        # never reassigned mat... so this sould be because of unknown extension.
        # shouldn't be able to reach this point.  ha ha ha
        raise RuntimeError('In <function read_matrix_from_file>, you reached'
                           'an unreachable point.  Congrats.')
    else:
        if return_array:
            if bln_got_df:
                # we have ids, in the form of the index and columns
                return mat.to_numpy()  #this drops the index and columns
            else:
                if matrix_has_ids:
                    # remove first row and column
                    return mat[1:, 1:]
                else:
                    return mat
        else:
            if matrix_has_ids:
                if bln_got_df:
                    # already got the df_ffs
                    pass
                else:
                    row_ids = try_convert_ids_to_int(mat[1:, 0])
                    col_ids = try_convert_ids_to_int(mat[0, 1:])

                    df = pd.DataFrame(data=mat[1:, 1:],
                                    index=row_ids, columns=col_ids)
            else:
                # get geo IDs to make a dataframe

                theIDs, _ = getput_idsnames_fromfiles(filenameID=geofile,
                                                      table=cfg.tbl_files,
                                                      return_dict=False)
                # create dataframe
                # note that we may create sparse matrices...
                # to convert: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sparse.to_dense.html
                # there are probably more elegant/faster/more efficient ways to do this.

                if isinstance(mat, scipy.sparse.csc.csc_matrix):
                    print(f'\t\tconverting sparse matrix {file + extension} to dataframe...')
                    df = pd.DataFrame.sparse.from_spmatrix(data=mat,
                                                           index=theIDs.array,
                                                           columns=theIDs.array)
                else:
                    print(f'\t\tconverting matrix {file + extension} to dataframe...')
                    df = pd.DataFrame(data=mat, index=theIDs.array, columns=theIDs.array)
                print('\t\t\t...done')

            return df



def getput_idsnames_fromfiles(filenameID, table, return_dict=False):
    """
    Pull out ids and names from the file table, and write ids and names into
    appropriate columns in the file tables
    We store ids, names as dictionaries...(to avoid errors putting dataframes or series into dataframes)
    # note we must use ordered dict:
    # https://stackoverflow.com/questions/18996714/how-to-turn-pandas-dataframe-row-into-ordereddict-fast

    # but there is an option to return either the series/dataframe or dict
    :param filenameID: row index in one of the data tables
    :param table: name of the table
    :param return_dict: boolean; return to calling function a dataframe, unless this is true
    :return: ids as series, names as dataframe (index is the ids)
    """

    # Get the associated geofile.
    geoFileNameID = ''

    s_no_name = 'No names for geofile'

    if table == cfg.tbl_files:
        if cfg.df_files.loc[filenameID, cfg.tbl_files_col_struct] == cfg.str_fileStruct_geo:
            geoFileNameID = filenameID
        else:
            geoFileNameID = cfg.df_files.loc[filenameID, cfg.tbl_files_col_assocGeo]

    elif table == cfg.tbl_data:
        geoFileNameID = cfg.df_data.loc[filenameID, cfg.tbl_data_col_assocGeo]

    else:
        raise KeyError(f'function <getput_idsnames_fromfiles> was passed the wrong kind of table')

    if geoFileNameID == '':
        raise KeyError(f'function <getput_idsnames_fromfiles> did not find filename = <{filenameID} in table = <{table}>')

    # now we know the geofile
    # first we check if names and IDs are already here
    datarow = cfg.df_files.loc[filenameID, :]

    # check if we have the IDs (names will be included, if available)
    bln_have_ids = filenameID in cfg.dict_ids.keys()

    if bln_have_ids:
        # we already have gotten the ids; can simply get them out of the dictionary
        pass
    else:
        # missing, so we need to get the file
        filename, path, ext = get_filepathext_from_files(filenameID)

        if not os.path.isfile(path + filename + ext):
            raise FileExistsError(f'Function <getput_idsnames_fromfiles> did not find file {path + filename + ext}, '
                                  f'and {filenameID} in table {table}')

        col_id = datarow[cfg.tbl_files_col_headerID]
        col_name = datarow[cfg.tbl_files_col_headerName]
        col_text2 = datarow[cfg.tbl_files_col_headerText2]
        bln_have_name_header = checkstr(col_name)
        bln_have_text_header2 = checkstr(col_text2)

        # return IDs as series; we don't need to fix_ids because we do it in this function.
        series_ids = read_column_from_file(path=path, file=filename, extension=ext,
                                           col=col_id,
                                           fix_IDs=False,
                                           return_array=False,
                                           realIDcolumn= col_id)

        # series_ids = read_column_from_file(path=path, file=filename, extension=ext,
        #                                    col=col_id, colID=col_id,
        #                                    fix_IDs=False,
        #                                    return_array=False)

        try:
            series_ids =series_ids.astype('int')
        except ValueError:
            pass

        #rename the series
        series_ids.rename(cfg.s_authoritative_ids, inplace=True)

        # now series_ids remains as it was read, without sorting, dropping duplicates, etc.

        if not bln_have_name_header:
            # create unique
            series_ids_unique = series_ids.drop_duplicates()
            series_ids_unique.sort_values(axis=0, ascending=True, inplace=True)
            series_ids_unique.reset_index(drop=True, inplace=True)  # drop to prevent the index from being added as a column, making it a dataframe

            # add both to the respective dictionaries
            cfg.dict_ids[filenameID] = series_ids
            cfg.dict_ids_unique[filenameID] = series_ids_unique

        if bln_have_name_header:
            # there is a name field, but we don't have names yet
            # we'll define the unique IDs and names based on the combination of the two.
            #   Note the IDs and names *should* be the same duplicates, but need to build in a check

            # don't need to fix IDs, because we do later in this function
            series_names = read_column_from_file(path=path, file=filename, extension=ext,
                                                 col=col_name,
                                                 fix_IDs=False,
                                                 return_array=False,
                                                 realIDcolumn= datarow[cfg.tbl_files_col_headerID])
            series_names.rename (cfg.s_authoritative_names, inplace = True)

            if bln_have_text_header2:
                series_text2 = read_column_from_file(path=path, file=filename, extension=ext,
                                                     col=col_text2,
                                                     fix_IDs=False,
                                                     return_array=False,
                                                     realIDcolumn=datarow[cfg.tbl_files_col_headerID])
                series_text2.rename(cfg.s_authoritative_names, inplace=True)

            # create dataframe with two columns, IDs and names.  Note that the index is not the IDs
            df_names = pd.DataFrame({cfg.s_authoritative_ids: series_ids, cfg.s_authoritative_names: series_names.values})
            if bln_have_text_header2:
                df_names[cfg.s_authoritative_ISO3] = series_text2
            # df_names = pd.Dataframe(index=series_ids, data=series_names)

            # the duplicates should be the same whether dropping duplicates in just IDs,  Names, or in both.
            #   I.e., we assume that each ID corresponds to one name.
            #   In case that's not true, we prioritize the IDs, so use *only* the IDs to drop duplicates

            df_names_unique = df_names.drop_duplicates(subset=[cfg.s_authoritative_ids])
            df_names_unique.sort_values(axis=0, by=cfg.s_authoritative_ids, ascending=True, inplace=True)
            df_names_unique.reset_index(drop=True, inplace=True)

            series_ids = df_names.loc[:,cfg.s_authoritative_ids]
            series_ids_unique = df_names_unique.loc[:,cfg.s_authoritative_ids]

            df_names = df_names.set_index(cfg.s_authoritative_ids)
            df_names_unique = df_names_unique.set_index(cfg.s_authoritative_ids)

            # add to the respective dictionaries.
            cfg.dict_ids[filenameID] = series_ids
            cfg.dict_ids_unique[filenameID] = series_ids_unique

            cfg.dict_names[filenameID] = df_names
            cfg.dict_names_unique[filenameID] = df_names_unique

    # now we have either seen the ids and names are in the df_files, or we have put them there
    # so get the variables out of the dictionary:
    ids = cfg.dict_ids_unique[filenameID]

    bln_have_names = filenameID in cfg.dict_names_unique.keys()
    if bln_have_names:
        names = cfg.dict_names_unique[filenameID]  #this will be a dataframe, with index as the ids
    else:
        names = s_no_name

    if return_dict:
        # I don't think this is used...
        pass
    else:
        return ids, names


    #         # we should only get names and IDs from geofiles...
    # if table == cfg.t_geo_files:
    #     # we assume that ONLY geofiles have the full list of names
    #     # so we set a boolean to allow recording of ids
    #     bln_record_ids = True
    # else:
    #     # we will not record the ids into the geofiles table, but will just return what we read
    #     bln_record_ids = False
    #
    # if filenameID in cfg.dict_excelTables[table].index:
    #
    #     # first we check if names and IDs are already here
    #     datarow = cfg.dict_excelTables[table].loc[filenameID, :]
    #
    #     if bln_record_ids:
    #         if all(x in datarow.index for x in [cfg.add_geo_id, cfg.add_geo_name]):
    #             # we already have the id and name columns, so there may be data
    #             pass
    #         else:
    #             # there are no columns with the id and names, so we add them
    #             cfg.dict_excelTables[table][cfg.add_geo_id] = ''
    #             cfg.dict_excelTables[table][cfg.add_geo_name] = ''
    #
    #             # re-retrieve the data row -- now it has the added columns?
    #             datarow = cfg.dict_excelTables[table].loc[filenameID, :]
    #
    #         # now, we know that the data table contains ids and names
    #         # check if there are lists, dictionaries, etc... but NOT empty strings
    #         bln_have_ids = datarow[cfg.add_geo_id] != ''
    #         bln_have_names = datarow[cfg.add_geo_name] != ''
    #
    #     else:
    #         # bln_record_ids is false, so we get data anew
    #         bln_have_ids = False
    #         bln_have_names = False
    #
    #     # now we have determined whether we already have names/ids in the geofile table
    #     # (Note: in the case of looking up a file in the regular file table,
    #     #   we DO NOT get ids from the geofile table, but rather read from the new file)
    #
    #     if not bln_have_ids or not bln_have_names:
    #         # one is missing, so we need to get the file
    #         filename, path, ext = get_filepathext_from_files(filenameID, table)
    #
    #         if not os.path.isfile(path + filename + ext):
    #             raise FileExistsError(f'Function oldget_data_from_datanameID '
    #                                   f'did not find file {path + filename + ext}, '
    #                                   f'and {filenameID} in table {table}')
    #
    #         col_id = datarow[cfg.xlcols[table][cfg.s_id]]
    #
    #         # for the geofiles table, there will be a name column, but not necessarily a value
    #         # for regular files, there is NOT a name column
    #         if table == cfg.t_geo_files:
    #             # this is a bit kludgy... this function was originally for just the
    #             # table geofiles, but I expanded it to include both tables.
    #              if cfg.s_geofile_namecol in datarow.index:
    #                  col_name = datarow[cfg.s_geofile_namecol]
    #         else:
    #             col_name = None
    #
    #         if col_name is None:
    #             bln_have_name_header = False
    #         else:
    #             if len(col_name) > 0:
    #                 bln_have_name_header = True
    #             else:
    #                 bln_have_name_header = False
    #
    #         # we'll need the ids whether or not we're missing ids or names
    #         # TODO: ids and names could return numbers represented as strings...
    #         # e.g., the GeosGrid from Roy.  Do we force conversion to number?
    #
    #         if table == cfg.t_reg_files:
    #             # need to set a sheet name
    #             xlsheet = datarow[cfg.s_regfile_xlsheet]
    #         else:
    #             xlsheet = ''
    #
    #         ids = read_column_from_file(path=path, file=filename, extension=ext,
    #                                     col=col_id, sheet=xlsheet,
    #                                     return_array=False)
    #         try:
    #             ids = ids.astype(int)
    #         except ValueError:
    #             # if we cannot convert to int, we get a value error
    #             pass
    #
    #         if bln_have_name_header and not bln_have_names:
    #             # there is a name field, but we don't have names yet
    #             names = read_column_from_file(path=path, file=filename,
    #                                           extension=ext, col=col_name)
    #
    #             # join the ids (as index) to the name col_values
    #             df_names = pd.DataFrame(data=names.values, index=ids.values)
    #
    #         if bln_record_ids:
    #             # write to cfg.dict_excelTables
    #
    #             # write the ids if we don't have them
    #             if not bln_have_ids:
    #                 cfg.dict_excelTables[table].loc[filenameID, cfg.add_geo_id] = \
    #                     [ids.to_dict(into=OrderedDict)]
    #
    #
    #             if bln_have_name_header:
    #                 if not bln_have_names:
    #                     cfg.dict_excelTables[table].loc[filenameID, cfg.add_geo_name] = \
    #                         [df_names.to_dict(into=OrderedDict)]
    #                 else:
    #                     pass
    #                     # names are already there, and we didn't read them, anyhow.
    #             else:
    #                 cfg.dict_excelTables[table].loc[filenameID,
    #                                         cfg.add_geo_name] = 'No Name'
    #
    #     # now we either had the IDS and names, or we have created them...
    #     # These are dictionaries... though the name might be 'no name'
    #     # so we can return:
    #     if not bln_record_ids:
    #         return ids
    #     else:
    #
    #         # get them back out of the table
    #         ids = cfg.dict_excelTables[table].loc[filenameID, cfg.add_geo_id]
    #         names = cfg.dict_excelTables[table].loc[filenameID, cfg.add_geo_name]
    #
    #         if return_dict:
    #             return ids, names
    #         else:
    #             if names == 'No Name':
    #                 return pd.Series(ids), names
    #             else:
    #                 # we have both ids and names that can be converted to series and dict
    #                 return pd.Series(ids), pd.DataFrame.from_dict(names,
    #                                                               orient='columns')


def getput_idsnames(datanameID, return_dict=False):
    # figure out which geofilenameID and table we're working in
    # geofilenameID is the unique index value in one of the
    # two data tables

    # returns ids and names (either as dfs or dicts)

    filenameID, table, _ = get_filetablecol(datanameID)

    # we can return the output of the file names directly
    return getput_idsnames_fromfiles(filenameID=filenameID,
                                     table=table,
                                     return_dict=return_dict)


