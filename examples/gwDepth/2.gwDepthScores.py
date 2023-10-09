from pcrglobwb_eval import grdc

#Directory containing the processed gw data (using 1.gwDepthPre.py)
gwDepthDataDirectory='/scratch/depfg/otoo0001/transient_analysis/gw_part/'

#Directory containing pcrglobwb output
simDirectory = '/scratch/depfg/otoo0001/transient_analysis/Nicole_dataset/'

gwDepthData = grdc(gwDepthDataDirectory, simDirectory)
validationData = gwDepthData.validationData
scores = gwDepthData.scores

print(validationData)
print(scores)