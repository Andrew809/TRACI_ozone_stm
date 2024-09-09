

from datetime import datetime

import calc.setup.config as cfg

from calc.setup.read_inputs import get_all_excel_tables
from calc.inout.files import make_filepathext, get_filepathext_from_files
from calc.setup.create_stmats import create_all_sdmats

from calc.functions.checks import checkstr

import pandas as pd

# display options for pandas output to console window
# pd.set_option('max_columns', 8)
# pd.set_option('max_colwidth', 40)
# pd.set_option('display.width', 200)

# some little checking functions


def main_setup(main_excel_file,
               stmat_file='', manual_stmat_skip=False,
               fftot_file='', save_fftots=False, manual_fftot_skip=False,
               save_fftot_native_shapefiles=False, save_fftot_calc_shapefiles=False,
               create_intersects=True,
               bln_normalize=False):

    # a global variable to append to outputs
    cfg.suffix_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

    # read all excel tables, populate the cfg.dict_excelTables dictionary
    cfg.dict_excelTables = get_all_excel_tables(fullfile= make_filepathext(file=main_excel_file))

    # assign the tables; # so to access a table, need to use cfg.dict_excelTables[tbl_data]
    cfg.df_files = cfg.dict_excelTables[cfg.tbl_files]
    cfg.df_data = cfg.dict_excelTables[cfg.tbl_data]
    cfg.df_stms = cfg.dict_excelTables[cfg.tbl_stms]
    cfg.df_flows = cfg.dict_excelTables[cfg.tbl_flows]
    cfg.df_impactScale = cfg.dict_excelTables[cfg.tbl_impactScale]
    cfg.df_calc = cfg.dict_excelTables[cfg.tbl_calc]

    # make tweaks to the tables

    # rename or copy columns we can use directly in output
    cfg.df_flows[cfg.str_outputCol_FlowName] = cfg.df_flows[cfg.tbl_flows_col_dispNameFull]
    cfg.df_flows[cfg.str_outputCol_FLowCAS] = cfg.df_flows[cfg.tbl_flows_col_CAStext]

    cfg.df_impactScale[cfg.str_outputCol_unit] = cfg.df_impactScale[cfg.tbl_impactScale_col_Unit]

    # files
    # empty columns that we will put IDs and Names into when we 1st read the files.
    # cfg.df_files[cfg.tbl_files_colADD_geoIDs] = ''
    # cfg.df_files[cfg.tbl_files_colADD_geoNames] = ''
    cfg.df_files[cfg.tbl_files_colAdd_pathFileExt] = ''

    # create a column of the full file location (path + file + ext)
    #   this comes in handy for looking up specific files from functions that have the full file location but not the unique ID
    for idx, row in cfg.df_files.iterrows():
        #print(calc_idx)
        if checkstr(row[cfg.tbl_files_col_fileName]):
            file, path, ext = get_filepathext_from_files(idx)
            cfg.df_files.at[idx, cfg.tbl_files_colAdd_pathFileExt] = (
                make_filepathext(file=file, path=path, ext=ext)
            )



    if manual_stmat_skip:
        pass
    else:
        if stmat_file != '':
            # this function takes only a file name
            create_all_sdmats(sdmat_store_file=stmat_file)
            # otherwise, sdmats are created on the fly



