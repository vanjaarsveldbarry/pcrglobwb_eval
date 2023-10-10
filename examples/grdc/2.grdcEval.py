from pcrglobwb_eval import grdc

#Directory containing the processed GRDC data (using 1.grdcPre.py)
grdcDataDirectory='/scratch/depfg/7006713/temp/validationDev/grdc'

#Directory containing pcrglobwb output
simDirectory = '/scratch/depfg/7006713/valDataTest/clone_13'

grdcData = grdc(grdcDataDirectory, simDirectory)
validationData = grdcData.validationData
scores = grdcData.scores

print(validationData)
print(scores)
