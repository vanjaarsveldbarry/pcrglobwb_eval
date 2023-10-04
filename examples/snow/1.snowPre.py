from pcrglobwb_eval import preprocess

#Directory containing the raw snow cover data:
snowRawDataDirectory='/scratch/depfg/7006713/data/validationData/snowData/RAW'

#Directory where to save the preprocessed snow cover data:
saveFolder='/scratch/depfg/7006713/temp/validationDev/snow'
preprocess.snowData(snowRawDataDirectory=snowRawDataDirectory, saveFolder=saveFolder)