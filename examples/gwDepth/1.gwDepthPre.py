from pcrglobwb_eval import preprocess

#Directory containing the GRDC raw data files(i.e. 1104150_Q_Day.Cmd.txt)
groundWaterDepthDataDirectory='/scratch/depfg/otoo0001/data/Groundwater_observation_data_Australia/AUS - Well and Monitoring Data/monitoring'


#Directory where to save the preprocessed GRDC data:
saveFolder='/scratch/depfg/otoo0001/transient_analysis/gw_part/'
# - GRDC_array.zarr
# - GRDC_points.shp

preprocess.groundWaterDepthData(groundWaterDepthDataDirectory=groundWaterDepthDataDirectory, saveFolder=saveFolder)