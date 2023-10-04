from pcrglobwb_eval import preprocess

#Directory containing the GRDC raw data files(i.e. 1104150_Q_Day.Cmd.txt)
grdcRawDataDirectory='/scratch/depfg/7006713/data/validationData/GRDC/_GRDC_raw/data'

#FilePath to High Resolution ldd:
lddFile='/scratch/depfg/7006713/temp/validationDev/grdc/lddsound_30sec_version_202005XX_correct_lat.nc'

#Directory where to save the preprocessed GRDC data:
saveFolder='/scratch/depfg/7006713/temp/validationDev/grdc/test'


preprocess.grdcData(grdcRawDataDirectory=grdcRawDataDirectory, 
                    saveFolder=saveFolder, 
                    lddFile=lddFile)