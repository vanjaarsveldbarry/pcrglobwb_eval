from pcrglobwb_eval import grdc

#Directory containing the processed GRDC data (using 1.grdcPre.py)
grdcDataDirectory='/scratch/depfg/7006713/temp/validationDev/grdc'

#Directory containing pcrglobwb output
simDirectory = '/scratch/depfg/7006713/temp/eguData/geowat_eu/30sec/30secNewDownscale/clone_3/1979_1989'

grdcData = grdc(grdcDataDirectory, simDirectory)
validationData = grdcData.validationData
scores = grdcData.scores

print(validationData)
print(scores)