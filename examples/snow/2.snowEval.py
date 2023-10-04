from pcrglobwb_eval import snow

#Directory containing the processed GRDC data (using 1.grdcPre.py)
snowDataDirectory='/scratch/depfg/7006713/temp/validationDev/snow'

#Directory containing pcrglobwb output
simDirectory = '/scratch/depfg/7006713/temp/eguData/geowat_eu/30sec/30secNewDownscale/clone_3/1990_2000'

snowData = snow(snowDataDirectory, simDirectory)
validationData = snowData.validationData
print(validationData)
scores = snowData.scores
print(scores)
# print(validationData)
# print(scores)