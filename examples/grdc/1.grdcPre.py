from pcrglobwb_eval import preprocess

#Directory containing the GRDC raw data files(i.e. 1104150_Q_Day.Cmd.txt)
grdcRawDataDirectory='/scratch/steya001/_GRDC_raw/data'

#Directory where to save the preprocessed GRDC data:
saveFolder='/scratch/steya001/test_eval'
preprocess.grdcData(grdcRawDataDirectory=grdcRawDataDirectory, saveFolder=saveFolder)
