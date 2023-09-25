from pcrglobwb_eval import grdc

#Directory containing the processed GRDC data
grdcDataDirectory='/scratch/depfg/7006713/data/validationData/GRDC'

#Directory containing pcrglobwb output
simDirectory = '/scratch/depfg/7006713/geowat_global_output/5min/5landcoverTypes/clone_01/1979_2019/'

grdcData = grdc(grdcDataDirectory, simDirectory)
validationData = grdcData.validationData
scores = grdcData.scores
print(validationData)
print(scores)