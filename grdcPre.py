from pcrglobwb_eval import preprocess

grdcRawDataDirectory='/scratch/depfg/7006713/data/validationData/GRDC/_GRDC_raw/data'
saveFolder='/scratch/depfg/7006713/temp/test'
preprocess.grdcData(grdcRawDataDirectory=grdcRawDataDirectory, saveFolder=saveFolder)