from pcrglobwb_eval import grdc

grdcObservedDataDirectory='/scratch/depfg/7006713/data/validationData/GRDC'
simDirectory = '/scratch/depfg/7006713/geowat_global_output/30sec/clone_1'
savePath = '/scratch/depfg/7006713/temp/test/output'

grdcData = grdc(grdcObservedDataDirectory, simDirectory)
validationData = grdcData.validationData
scores = grdcData.scores
