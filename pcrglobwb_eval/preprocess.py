
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import warnings
import numpy as np
import xarray as xr
import shutil
import geopandas as gpd
from functools import partial
import logging
import concurrent.futures
from shapely.geometry import Point
import time

def timing_decorator(func):
    import time
    '''@timing_decorator ontop of function you want to time'''
   
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"{func.__name__} took: {execution_time:.6f} seconds to run")
        return result
    


logging.basicConfig(filename="preprocessing.log", level=logging.INFO)

class preprocess:
    """Preprocessig input data for further validation.
    """
   
    def groundWaterDepthData(groundWaterDepthDataDirectory, saveFolder):

        saveFolder = Path(saveFolder)
        groundWaterDepthDataDirectory = Path(groundWaterDepthDataDirectory)

        saveFolder.mkdir(exist_ok=True)

        
        import concurrent.futures
        def makeArrayDatasetGW(stationFiles):
            (saveFolder / "temp").mkdir(exist_ok=True)
            stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx"))  
            totalLength = len(stationFiles)
            count = 0
            stationFiles = stationFiles
            dataSets = []
            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = []
                for file in tqdm(stationFiles):
                    futures.append(executor.submit(makeArrayDatasetGW, file))
                for future in concurrent.futures.as_completed(futures):
                    dataFile = future.result()
                    count = count + 1
                    selected_columns = ['ID', 'Date and Time', 'Value']
                    dataFile = pd.read_excel(file, skiprows=[1], usecols=selected_columns)
                    dataFile['Date and Time'] = pd.to_datetime(dataFile['Date and Time']).dt.date
                    dataFile = dataFile.dropna()
                    if dataFile.size <= 10: continue

                    if count == 1:
                        wellFilePath=file.parents[1] / 'wells.xlsx'
                        welldataFile = pd.read_excel(wellFilePath, skiprows=[1], usecols=['ID', 'Latitude', 'Longitude'])

                    wellID = dataFile['ID'].iloc[0]
                    locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]

                    if locationData.empty:
                        print(f'Well {wellID} is not found') 
                        continue

                    data = dataFile["Value"]
                    data = data.values.reshape(data.shape + (1,))

                    timeVals = pd.to_datetime(dataFile["Date and Time"])
                    wellVals = np.array([wellID]).astype("str")

                    dataSet = xr.Dataset(
                        data_vars={
                            "gwDepth": ((["time", "well"]), data),
                        },
                        coords={"time": (["time"], timeVals),
                                "well": (["well"], wellVals),
                                "lat": (["well"], np.array([locationData["Latitude"].iloc[0]]).astype("float32")),
                                "lon": (["well"], np.array([locationData["Longitude"].iloc[0]]).astype("float32"))})

                    dataSet["gwDepth"] = dataSet["gwDepth"].dropna(dim="time")
                    dataSet = dataSet.drop_duplicates(dim='time')

                    dataSets.append(dataSet)
                
            full_ds = xr.concat(dataSets, dim="well")
            full_ds = full_ds.sortby("lat", ascending=False)
            full_ds.to_zarr(saveFolder / "gwDepth_array.zarr", mode="w", consolidated=True)
            logging.info(f"zaar saved sucessfully: {file}")
            
        def makedataPointsfromzarr(saveFolder):
            zarr_file = Path(saveFolder) / "gwDepth_array.zarr"
            zarr_data = xr.open_zarr(zarr_file)
            latitudes = zarr_data['lat'].values
            longitudes = zarr_data['lon'].values
            field_values = zarr_data["well"].values

            geometry = [Point(lon, lat) for lon, lat in zip(longitudes, latitudes)]
            (saveFolder / "well_points").mkdir(exist_ok=True)
            gdf = gpd.GeoDataFrame({"well": field_values, 'geometry': geometry}, crs="EPSG:4326")
            gdf.to_file(saveFolder / "well_points/obswell_points.shp")


        saveFolder = Path(saveFolder)
        groundWaterDepthDataDirectory= Path(groundWaterDepthDataDirectory)
        stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx")) 
        makeArrayDatasetGW(stationFiles)
        makedataPointsfromzarr(saveFolder)

        











                    # @timing_decorator
                    # def write_chunk_file(file, dataSet):
                    #     full_ds = dataSet(file)
                    #     stationIndex = np.where(dataSet.well.values == full_ds.well)[0][0]
                    #     timeStartIndex = np.where(dataSet.time.values == full_ds['ts_start'])[0][0]
                    #     timeEndIndex = np.where(dataSet.time.values == full_ds['ts_end'])[0][0]
                    #     full_ds.to_zarr(saveFolder / "gwDepth_array.zarr", region={"time": slice(timeStartIndex, timeEndIndex+1), "well": slice(stationIndex,stationIndex+1)})
                    #     print(full_ds)

                    # with concurrent.futures.ProcessPoolExecutor(max_workers=80) as executor:
                    #     futures = []
                    #     for file in stationFiles:
                    #         futures.append(executor.submit(write_chunk_file, file, dataSet))
                    #     for future in tqdm(concurrent.futures.as_completed(futures), total=len(stationFiles), desc="Writing zarr files"): pass
           



          # zarr_file = Path(saveFolder) / "gwDepth_array.zarr"
            # zarr_data = xr.open_zarr(zarr_file)
            # latitudes = zarr_data['lat'].values
            # longitudes = zarr_data['lon'].values
            # field_values = zarr_data["well"].values

            # geometry = [Point(lon, lat) for lon, lat in zip(longitudes, latitudes)]
            # (saveFolder / "well_points").mkdir(exist_ok=True)
            # gdf = gpd.GeoDataFrame({"well": field_values, 'geometry': geometry}, crs="EPSG:4326")
            # gdf.to_file(saveFolder / "well_points/obswell_points.shp")














            # with concurrent.futures.ProcessPoolExecutor(max_workers=80) as executor:
            #     futures = []
            #     for file in stationFiles:
            #         futures.append(executor.submit(makeArrayDatasetGW, file))
            #         for future in tqdm(concurrent.futures.as_completed(futures)): pass
            #     full_ds = futures[0].result()

            # print(full_ds)
 
            # def save_zarr(full_ds):
        #     full_ds.to_zarr(saveFolder / "gwDepth_array.zarr", mode="w", consolidated=True)
        #     logging.info(f"zaar saved sucessfully: {file}") 


# import concurrent.futures


# def timing_decorator(func):
#     import time
#     '''@timing_decorator ontop of function you want to time'''

#     def wrapper(*args, **kwargs):
#         start_time = time.time()
#         result = func(*args, **kwargs)
#         end_time = time.time()
#         execution_time = end_time - start_time
#         print(f"{func.__name__} took: {execution_time:.6f} seconds to run")
#         return result
    
    
#     return wrapper
# logging.basicConfig(filename="preprocessing.log", level=logging.INFO)
# class preprocess:
#     """Preprocessig input data for further validation.
#     """
#     def groundWaterDepthData(groundWaterDepthDataDirectory, saveFolder):

#         saveFolder = Path(saveFolder)
#         groundWaterDepthDataDirectory = Path(groundWaterDepthDataDirectory)

#         saveFolder.mkdir(exist_ok=True)

#         def makeArrayDatasetGW(file):
#             selected_columns = ['ID', 'Date and Time', 'Value']
#             dataFile = pd.read_excel(file, skiprows=[1], usecols=selected_columns)
#             dataFile['Date and Time'] = pd.to_datetime(dataFile['Date and Time']).dt.date
#             dataFile = dataFile.dropna()
#             if dataFile.size <= 100: return None
            
#             wellFilePath=file.parents[1] / 'wells.xlsx'
#             welldataFile = pd.read_excel(wellFilePath, skiprows=[1], usecols=['ID', 'Latitude', 'Longitude'])

#             wellID = dataFile['ID'].iloc[0]
#             locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]

#             if locationData.empty:
#                 print(f'Well {wellID} is not found') 
#                 return None

#             data = dataFile["Value"]
#             data = data.values.reshape(data.shape + (1,))

#             timeVals = pd.to_datetime(dataFile["Date and Time"])
#             wellVals = np.array([wellID]).astype("str")

#             dataSet = xr.Dataset(
#                 data_vars={
#                     "gwDepth": ((["time","well"]), data),
#                 },
#                 coords={"time": (["time"], timeVals),
#                         "well": (["well"], wellVals),
#                         "lat": (["well"], np.array([locationData["Latitude"].iloc[0]]).astype("float32")),
#                         "lon": (["well"], np.array([locationData["Longitude"].iloc[0]]).astype("float32"))})

#             dataSet["gwDepth"] = dataSet["gwDepth"].dropna(dim="time")
#             dataSet = dataSet.drop_duplicates(dim='time')

#             return dataSet

#         def save_zarr(full_ds):
#             full_ds.to_zarr(saveFolder / "gwDepth_array.zarr", mode="w", consolidated=True)
#             logging.info(f"zaar saved sucessfully: {file}")

#         stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx"))[:10000]
#         with concurrent.futures.ProcessPoolExecutor() as executor:
#             futures = []
#             for file in stationFiles:
#                 futures.append(executor.submit(makeArrayDatasetGW, file))
#             full_ds = None
#             for future in concurrent.futures.as_completed(futures):
#                 result = future.result()
#                 if result is not None:
#                     if full_ds is None:
#                         full_ds = result
#                     else:
#                         full_ds = xr.merge([result, full_ds], join="outer", compat="no_conflicts")
#                         full_ds = full_ds.sortby("lat", ascending=False)
#         print(full_ds)
#         save_zarr(full_ds)

#         def makedataPointsfromzarr(saveFolder):
#             zarr_file = Path(saveFolder) / "gwDepth_array.zarr"
#             zarr_data = xr.open_zarr(zarr_file)
#             latitudes = zarr_data['lat'].values
#             longitudes = zarr_data['lon'].values
#             field_values = zarr_data["well"].values

#             geometry = [Point(lon, lat) for lon, lat in zip(longitudes, latitudes)]
#             (saveFolder / "well_points").mkdir(exist_ok=True)
#             gdf = gpd.GeoDataFrame({"well": field_values, 'geometry': geometry}, crs="EPSG:4326")
#             gdf.to_file(saveFolder / "well_points/obswell_points.shp")


#         saveFolder = Path(saveFolder)
#         groundWaterDepthDataDirectory= Path(groundWaterDepthDataDirectory)
#         stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx")) 
#         makeArrayDatasetGW(stationFiles)
#         makedataPointsfromzarr(saveFolder)









 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#         # try:
#         #     saveFolder = Path(saveFolder)
#         #     groundWaterDepthDataDirectory = Path(groundWaterDepthDataDirectory)

#         #     saveFolder.mkdir(exist_ok=True)
#         #     count = 0
#         #     stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx"))
#         #     full_ds = None  # Initialize the full dataset

#         #     for file in tqdm(stationFiles):
#         #         count = count + 1
#         #         # Load station info from the 'wells.xlsx' file
#         #         wellFilePath = file.parents[1] / 'wells.xlsx'
#         #         welldataFile = preprocess.read_well_data(wellFilePath)
#         #         welldataFile = welldataFile[:1]

#         #         if welldataFile is None:
#         #             continue
#         #         # Extract station info from 'welldataFile'
#         #         stationInfo = welldataFile  

#         #         stationgGrid = np.array([stationInfo["ID"]])
#         #         stationgGrid = stationgGrid.reshape((1,) + stationgGrid.shape)

#         #         dataSet = pd.DataFrame(
#         #             {
#         #                 "stationID": np.array([stationInfo["ID"]]),
#         #                 "lat": np.array([stationInfo["Latitude"]]),
#         #                 "lon": np.array([stationInfo["Longitude"]]),
#         #             }
#         #         )
#         #         if full_ds is None:
#         #             full_ds = dataSet
#         #         else:
#         #             full_ds = pd.concat([dataSet, full_ds], axis=0, ignore_index=False)

#         #     (saveFolder / "obswells").mkdir(exist_ok=True)
#         #     gdf = gpd.GeoDataFrame(full_ds, geometry=gpd.points_from_xy(full_ds.lon, full_ds.lat), crs="EPSG:4326")
#         #     gdf.to_file(saveFolder / "obs_wells/wells.shp")
#         #     logging.info("Shapefile created: %s", saveFolder / "obs_wells/wells.shp")
#         # except Exception as e:
#         #     error_message = f"Error in makePointsDataset: {str(e)}"
#         #     logging.error(error_message)
#         #     print(error_message)