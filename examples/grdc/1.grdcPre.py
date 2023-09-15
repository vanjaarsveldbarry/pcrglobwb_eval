from pcrglobwb_eval import preprocess

#Directory containing the GRDC raw data files(i.e. 1104150_Q_Day.Cmd.txt)
grdcRawDataDirectory='/scratch/depfg/7006713/data/validationData/GRDC/_GRDC_raw/data'

#Directory where to save the preprocessed GRDC data as \:
# - GRDC_array.zarr
# - GRDC_points.shp

saveFolder='/scratch/depfg/7006713/temp/test'

preprocess.grdcData(grdcRawDataDirectory=grdcRawDataDirectory, saveFolder=saveFolder)