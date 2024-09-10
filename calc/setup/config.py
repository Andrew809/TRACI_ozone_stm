# https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
# https://stackoverflow.com/questions/13034496/using-global-variables-between-files

import os

# region User Settings

#---- Aggregation with weights
# aggregation with weights: what to do about targets with NO valid weighting values?
#   Note that a missing value is a NaN, not a zero.  A zero gives an area a 0 weight.
#   this boolean is applied only in the case where an entire target has no weight.
#   E.g., a small country with no fertilizer application data.
bln_area_weight_if_needed = True

# ---- Sector weighting
# calculate general (unknown) sector weighting if possible
#   if both agricultural and non-agricultural are provided, calculate the general.
#   TODO: define other methods of this calculation
#   For now, we use an average of the two.
bln_calc_sector_general = True

# HARD-CODED: these must match excel, in Tbl.Calc, also in list.Weight_sector
str_sector_agric = 'Agricult'
str_sector_notAg = 'NotAgri'
str_sector_genrl = 'General'


# what do we report in the final output table?  
# (Note that in all (some?) cases, we need to calculate agricultural and non-agric. to calculate general)
dict_sector_report = {str_sector_agric:False, str_sector_notAg:False, str_sector_genrl:True}


# ------- GEP aggregation
# must match Excel - HARD-CODED
str_GEPsetup_normalize_across_receiving = 'normalize across receiving elements'
str_GEPsetup_normalize_entire_matrix = 'normalize the full [emit x receive] matrix'
str_GEPsetup_weight_rows_by_vector = 'weight receiving rows by vector'

str_GEPsetup_selected = str_GEPsetup_normalize_entire_matrix

# endregion

# a global suffix for saving files
suffix_time = ''  # this will be set by datetime upon execution

# debugging levels
bln_debug = True
bln_debug2 = False


# region File locations ----------------
#   where to store intermediate calculations,
#   e.g., the intersections
#   without a drive letter, paths are relative to package
#   NO trailing file separators...  this is HARD-CODED

# CHOOSE APPROPRIATE LOCATION HERE...
#dir_storage = r'\datasets\Intermed'
# # ADH home computer:
# dir_storage = r'C:\CalcDir\TRACI_ozone_stm\tempStorage'
# ADH office:
dir_storage = r'C:\CalcDir\TRACI_ozone_stm\tempStorage'
# # ML computer:
# dir_storage = r'C:\Users\mariole\OneDrive - NTNU\LC-IMPACT\Terrestrial_Acidification\CFs\ADH Code\Input Data'

dir_output = r'\datasets\Output'

dir_stms = dir_storage + r'\STMats'
dir_shapefiles = dir_storage + r'\Shapefiles'
dir_intersects = dir_storage + r'\Intersects'
bln_save_intersects = True  #big time-saver, as the multi-intersects can take a long time. (much longer in python than in arc or qgis)

dir_input_pkls = dir_storage + r'\PickledInputs'
bln_save_geopickles = True  #this speeds up reading the shapefiles a bit...

os.makedirs(dir_stms, exist_ok=True)
os.makedirs(dir_shapefiles, exist_ok=True)
os.makedirs(dir_intersects, exist_ok=True)
os.makedirs(dir_input_pkls, exist_ok=True)



str_fileText_Orig = 'Orig'
str_fileText_New = 'New'
str_fileText_Weight = 'Weight'

# endregion

# region ---- Output column names
str_outputCol_FlowName = 'FLOW_name'
str_outputCol_FLowCAS = 'FLOW_casnumber'
str_outputCol_lciaMethodLocISO3 = 'LCIAMethod_location'
str_outputCol_lciaMethodLocName = 'LCIAMethod_location_name'
str_outputCol_CF = 'CF'
str_outputCol_unit = 'Unit'
str_outputCol_CF_uncert_low = 'CF_Uncertainty_Lower'
str_outputCol_CF_uncert_high ='CF_Uncertainty_Higher'
str_outputCol_flowclass0 = 'FLOW_class0'
str_outputCol_flowclass1 = 'FLOW_class1'
str_outputCol_flowclass2 = 'FLOW_class2'
str_outputCol_species = 'Species'
str_outputCol_lciaMethodRealm = 'LCIAMethod_realm'
str_outputCol_lciaMethodMathAppr = 'LCIAMethod_mathematicalApproach'
str_outputCol_scenario = 'Scenario'
str_outputCol_CF_deriv = 'CF_derivation'
str_outputCol_matchCF = 'Matching_CF'
str_outputCol_matchComp = 'Matching_Compartment'
str_outputCol_emitResOrig = 'Emission Resolution Orig'
str_outputCol_emitResNew = 'Emission Resolution New'

list_out_cols = [
    str_outputCol_FlowName,
    str_outputCol_FLowCAS,
    str_outputCol_lciaMethodLocISO3,
    str_outputCol_lciaMethodLocName,
    str_outputCol_CF,
    str_outputCol_unit,
    str_outputCol_CF_uncert_low,
    str_outputCol_CF_uncert_high,
    str_outputCol_flowclass0,
    str_outputCol_flowclass1,
    str_outputCol_flowclass2,
    str_outputCol_species,
    str_outputCol_lciaMethodRealm,
    str_outputCol_lciaMethodMathAppr,
    str_outputCol_scenario,
    str_outputCol_CF_deriv,
    str_outputCol_matchCF,
    str_outputCol_matchComp,
    str_outputCol_emitResOrig,
    str_outputCol_emitResNew
    ]


# endregion

import pandas as pd
#pandas options: https://pandas.pydata.org/docs/user_guide/options.html

pd.options.mode.copy_on_write = True  #https://pandas.pydata.org/pandas-docs/stable/user_guide/copy_on_write.html#migrating-to-copy-on-write

# display
pd.options.display.max_columns = 15
pd.options.display.max_colwidth = 40
pd.options.display.width = 200
# pd.set_option('max_columns', 15)
# pd.set_option('max_colwidth', 40)
# pd.set_option('display.width', 200)








# region ---------- Shared dictionaries, dataframes:

# dictionary of tables from excel
dict_excelTables = {}

# dictionary of spatial distribution matrices
dict_STMats = {}

# dicitonary of RFs, EFs, and GEPs.
#   although there are <= 3 files for each, this'll save a bit of time in the execution, since we'll avoid re-reading files
dict_RFs = {}
dict_EFs = {}
dict_GEPs = {}
dict_GEPs_normalized = {}


# dictionaries of files we've already read -- note that geodata 'geometry' column is not kept, for space
#   keys are the same as the file unique name
bln_keep_filesHaveBeenRead = True
dict_filesHaveBeenRead = {}


# dictionaries of the authoritative series of IDs and Names for geofiles.
#   authoritative = we come back to these lists to make sure any data related to a geofile has the proper dimensions and order

s_authoritative_ids = 'IDs'
s_authoritative_names = 'Names'
s_authoritative_ISO3 = 'ISO3'

# dictionaries of values as read from the source
dict_ids = {} # will be series
dict_names = {}  #will dataframes, with col 1 = ids, and col 2 = names
# dictionaries with duplicate values removed
dict_ids_unique = {}
dict_names_unique = {}


# output multiindex (mi) names:
# for the columns
multiindex_col_emitResOrig = 'Original Emission Resolution'
multiindex_col_emitResNew = 'New Emission (Aggregation) Resolution'  # possible values set in tbl_calc
multiindex_col_subst = 'Substance'  # defined in tbl_calc
multiindex_col_impactScale = 'Impact Scale'  # defined it tbl_impScales
multiindex_col_sectorWeight = 'Sector Weight'  # defined in tbl_calc
multiindex_col_dataItem = 'Data Item' # possible values set below, using some of the output names
list_multiindex_col_levels = [multiindex_col_emitResOrig, multiindex_col_emitResNew, multiindex_col_subst, multiindex_col_impactScale, multiindex_col_sectorWeight, multiindex_col_dataItem]

# for the main index
multiindex_row_EmitIDs = 'Emit IDs'

multiindex_col_dataItem_CF = str_outputCol_CF
multiindex_col_dataItem_Name = str_outputCol_lciaMethodLocName
multiindex_col_dataItem_Code = str_outputCol_lciaMethodLocISO3

# create the blank dataframe with multiindex columns, and a single index for the rows
df_out_multiindex = pd.DataFrame(columns=pd.MultiIndex.from_tuples([], names=list_multiindex_col_levels),
                                 index=pd.Index([], name='Generic Index'))

dict_output_by_emitRes = {} # dictionary to hold the df_outs, indexed by emission resolution.

df_out_final = pd.DataFrame()  # declare this here, so it's available across modules

# endregion





# lots of names of tables, headings, and strings from excel.
# These are captured here so that if they are changed in excel, they need only be changed here.
# these could be stored in dictionaries, but then accessing them requires one more identifier.
# There are limited tables, so easiest to just access them with tbl_xxx

# in the case of geo files (shp, gpkg), the col_ID and col_names become the authoritative list of IDs, which are used to make sure that all matrix and vector operations keep data in the correct order

# Table of files
tbl_files = 'Tbl.Files'
df_files = pd.DataFrame()
tbl_files_col_uniqueName = 'File unique name'
tbl_files_col_struct = 'File structure'
tbl_files_col_fileDirMain = 'File directory main'
tbl_files_col_fileDirSub = 'File directory subfolder'
tbl_files_col_fileName = 'File name'
tbl_files_col_fileExt = 'File extension'
tbl_files_col_assocGeo = 'Associated geofile'
tbl_files_col_headerID = 'Header with ID'
tbl_files_col_headerName = 'Header with Name'
tbl_files_col_headerText2 = 'Header with Text 2'
#tbl_files_colADD_geoIDs = 'Geo IDs'
#tbl_files_colADD_geoNames = 'Geo Names'
tbl_files_colAdd_pathFileExt = 'Path + File + Extension'

# Table of data
# separated from files, since some data files could have multiple columns that we used
tbl_data = 'Tbl.Data'
df_data = pd.DataFrame()
tbl_data_col_uniqueName = 'Data unique name'
tbl_data_col_sourceUniqueName = 'Source file unique name'
tbl_data_col_excelSheetName = 'Excel sheet name'
    # assume that on this sheet, a table starts at A1
tbl_data_col_headerID = 'Header with ID'
tbl_data_col_headerData = 'Header with data'
tbl_data_col_units = 'Data units'
tbl_data_col_assocGeo = 'Associated geofile'
    # associated_geo: for non-geo files, this tells which authoritative IDs are associated with the data.  (For example, the data may not have values for all of the IDs in the authoritative list.)
tbl_data_col_dataType = 'Data type'
tbl_data_col_weightType = 'Weight type'
tbl_data_col_isMatrix = 'Is Matrix'
tbl_data_col_matHasIDs = 'Matrix has IDs'
    # matrix_has_ids: if true, we assume a 1st row and 1st column with IDs, such that data start at position 2,2 of the incoming data.  If false, we assume data starts at position 1,1.
tbl_data_col_matlabName = 'Matrix matlab name'
    # matrix_matlab_name: further specification for opening the matrix in matlab

# Table of spatial transforms
tbl_stms = 'Tbl.SpatialTransforms'
df_stms = pd.DataFrame()
tbl_stms_col_stmUniqueName = 'STM unique name'
tbl_stms_col_create = 'Create'
tbl_stms_col_geoOrig = 'Orig Geofile'
tbl_stms_col_geoNew = 'New Geofile'
tbl_stms_col_weightDataUnique = 'Weight Data'
tbl_stms_col_STMtype = 'STM type'
tbl_stms_col_fileDirMain = 'File directory main'
tbl_stms_col_fileDirSub = 'File directory subfolder'
tbl_stms_col_fileName = 'File name'
tbl_stms_col_fileExt = 'File extension'
tbl_stms_col_fileNameAuto = 'Auto file name'
tbl_stms_col_useExisting = 'Use Existing'

# Table of flows
tbl_flows = 'Tbl.Flows'
df_flows = pd.DataFrame()
tbl_flows_col_uniqueName ='Flow unique Name'
tbl_flows_col_runIt = 'Run'
tbl_flows_col_dispNameFull = 'Display name full'
tbl_flows_col_dispNameShort = 'Display name short'
tbl_flows_col_assocFF = 'Associated FF'
tbl_flows_col_CASnum = 'CAS num'
tbl_flows_col_CAStext = 'CAS text'

# Impact Scale Table
tbl_impactScale = 'Tbl.Impact_Scales'
df_impactScale = pd.DataFrame()
tbl_impactScale_col_Include = 'Include in Calc?'
tbl_impactScale_col_Unit = 'Unit'

# Match excel: HARD-CODED:
str_impactScale_regional = 'Regional'
str_impactScale_global = 'Global'


# Calculation table
tbl_calc = 'Tbl.Calc'
df_calc = pd.DataFrame()
tbl_calc_col_calcNum = 'Calc Num'
tbl_calc_col_doCalc = 'Do Calc?'
tbl_calc_col_substance = 'Substance'
tbl_calc_col_FF = 'FF'
tbl_calc_col_emitResOrig = 'Orig Emit Resolution'
tbl_calc_col_emitResNew = 'New Emit Resolution'
tbl_calc_col_emitWeightSector = 'Emit Weighting Sector'
tbl_calc_col_emitWeightData = 'Emit Weighting Data'
tbl_calc_col_emitSTM = 'Emit STM'
tbl_calc_col_recvResOrig = 'Orig Receive Resolution'
tbl_calc_col_recvResNew = 'New Receive Resolution'
tbl_calc_col_recvWeight = 'Receive Weighting'
tbl_calc_col_recvSTM = 'Receive STM'
tbl_calc_col_RF = 'RF'
tbl_calc_col_EF = 'EF'
tbl_calc_col_GEP = 'GEP'
tbl_calc_col_GEPnorm_method = 'GEP Normalization'


# these strings must match excel.  See source excel file, named ranges "list.x"
# HARD-CODED / MATCH
# some strings

# file structure types
str_fileStruct_geo = 'GIS (spatial)'
str_fileStruct_table = 'Table'
str_fileStruct_matrix = 'Matrix'

# data type
str_dataType_FF = 'Fate factor'
str_dataType_RF = 'Response factor'
str_dataType_EF = 'Effect factor'
str_dataType_GEP = 'GEP'
str_dataType_EFGEP = 'EF_GEP'
str_dataType_geom = 'Geometry'
str_dataType_weight = 'Weighting'

#weighting type
str_weightTypeinten = 'Intensive'
str_weightTypeexten = 'Extensive'

#STM types  (again, these must match excel... HARD-CODED)
str_STMtype_receive = 'newReceive'
str_STMtype_emit = 'newEmit'




# GEO
# names used for dealing with projections
proj_crs_wkt = 'well-known text'
proj_crs_code = 'authority code'
proj_crs_proj = 'proj'

proj_crs_default = proj_crs_wkt
proj_s_default = r'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["standard_parallel_1",0.0],UNIT["Meter",1.0]]'
proj_conv_default = 1 / 1e6  # meters, so we multiply by this to get km2