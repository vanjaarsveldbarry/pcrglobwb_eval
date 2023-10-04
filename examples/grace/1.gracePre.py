from pcrglobwb_eval import preprocess

#Directory containing the raw snow cover data:
graceRawDataDirectory='/scratch/depfg/7006713/data/validationData/GRACE_data/grace_downloaded_2023_05_may'

#Directory where to save the preprocessed snow cover data:
saveFolder='/scratch/depfg/7006713/temp/validationDev/grace'
preprocess.snowData(graceRawDataDirectory=graceRawDataDirectory, saveFolder=saveFolder)