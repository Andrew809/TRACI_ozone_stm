# set up calculation inputs

# define functions to read from the various excel tables
# Note that all column names are set in config.py

# read in all the tables (we will not use all)
# set up dictionaries to hold the tables and their headers.

import os  # check for file existence, use os.sep
from collections import OrderedDict

import geopandas as gpd
import numpy as np

import openpyxl  # table functionality requires v3.0.4 (or higher)


import pandas as pd
import scipy
from scipy import io


from calc.inout.files import add_full_directory


def checkstr(x): return x is not None and len(x) > 0


def read_excel_named_table(sheet, tablename,
                           firstColIntoIndex=True, firstRowIntoColumns=True):
    """
    Read a named table (different than a named range) from excel, and return as data frame
    :param sheet: an openpyxl worksheet object
    :param tablename: a string that is the table name
    :param firstColIntoIndex: boolean,
    do we convert the first column into the dataframe index?
    :param firstRowIntoColumns: boolean,
    do we convert the first row into the datafrom columns?
    :return: pandas dataframe; any row with all nas is dropped

    note: requires openpyxl v. 3.0.4 or higher
    """
    region = sheet.tables[tablename].ref

    if firstRowIntoColumns:
        df = pd.DataFrame(([cell.value for cell in row] for row in sheet[region][1:]),
                          columns=([cell.value for cell in sheet[region][0]]))
    else:
        df = pd.DataFrame(([cell.value for cell in row] for row in sheet[region][1:]))

    if firstColIntoIndex:
        df.set_index(df.columns[0], inplace=True)

    df.dropna(axis=0, how='all', inplace=True)

    return df


def get_all_excel_tables(fullfile):
    """
    Read excel CF file and get all the tables, return them as a dictionary
    :param fullfile: path to excel file with CalcCF info.
    :return: dictionary of tables; if there is a replacement name defined,
    this will be the key; otherwise, the excel table name is the key
    """

    wb = openpyxl.load_workbook(filename=fullfile, data_only=True)
    # print(wb.worksheets)

    dict_tables = {}

    for ws in wb.worksheets:
        for table in ws.tables.items():
            # table returns a tuple of name and address

            key_name = table[0]

            dict_tables[key_name] = read_excel_named_table(sheet=ws, tablename=table[0],
                                                           firstColIntoIndex=True,
                                                           firstRowIntoColumns=True)
    # print(dict_tables.keys())
    return dict_tables


