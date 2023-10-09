from tqdm import tqdm
from pathlib import Path
import pandas as pd
import warnings
import numpy as np
import xarray as xr
import shutil
import geopandas as gpd
import logging
from shapely.geometry import Point

logging.basicConfig(filename="preprocessing.log", level=logging.INFO)
class preprocess:
    """Preprocessig input data for further validation.
    """
    def groundWaterDepthData(groundWaterDepthDataDirectory, saveFolder):

        saveFolder = Path(saveFolder)
        groundWaterDepthDataDirectory = Path(groundWaterDepthDataDirectory)

        saveFolder.mkdir(exist_ok=True)

        def makeArrayDatasetGW(station_Files):
            (saveFolder / "temp").mkdir(exist_ok=True)
            stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx"))  
            totalLength = len(stationFiles)
            count = 0
            stationFiles = stationFiles
            for file in tqdm(stationFiles):
            # for file in stationFiles:
                
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
            
                
                # locationData_df = pd.DataFrame(locationData)
                # locationData_df.to_csv = (saveFolder / "location_data.csv")
                                               
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

                

                if count == 1:
                    full_ds = dataSet
                    

                else:
                    full_ds = xr.merge([dataSet, full_ds], join="outer", compat="no_conflicts")
                    full_ds = full_ds.sortby("lat", ascending=False)
        
            full_ds = full_ds.chunk({"well": 1})
            print(full_ds)
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
    
































# import logging
# from tqdm import tqdm
# from pathlib import Path
# import pandas as pd
# import numpy as np
# import xarray as xr
# import geopandas as gpd
# import matplotlib.pyplot as plt
# from shapely.geometry import Point

# # Configure logging
# logging.basicConfig(filename="preprocessing.log", level=logging.INFO)

# class preprocess:
#     """Preprocessing input data for further validation."""

#     @staticmethod
#     def read_excel_file(file):
#         """Read an Excel file and preprocess it."""
#         selected_columns = ["ID", "Date and Time", "Value"]
#         try:
#             dataFile = pd.read_excel(
#                 file, skiprows=[1], usecols=selected_columns, parse_dates=["Date and Time"]
#             )
#             dataFile["Date and Time"] = pd.to_datetime(
#                 dataFile["Date and Time"]
#             ).dt.date
#             if dataFile.shape[0] < 10:
#                 logging.warning(f"File {file} contains fewer than 10 rows of data. Skipping.")
#                 return None  # Skip this file
#             return dataFile
#         except Exception as e:
#             error_message = f"Error reading Excel file {file}: {str(e)}"
#             logging.error(error_message)
#             print(error_message)
#             return None

#     @staticmethod
#     def read_well_data(wellFilePath):
#         """Read well data from a separate file."""
#         try:
#             welldataFile = pd.read_excel(
#                 wellFilePath, skiprows=[1], usecols=["ID", "Latitude", "Longitude"]
#             )
#             return welldataFile
#         except Exception as e:
#             error_message = f"Error reading well data file {wellFilePath}: {str(e)}"
#             logging.error(error_message)
#             print(error_message)
#             return None
        
#     def groundWaterDepthData(groundWaterDepthDataDirectory, saveFolder):
#         saveFolder = Path(saveFolder)
#         groundWaterDepthDataDirectory = Path(groundWaterDepthDataDirectory)
#         saveFolder.mkdir(exist_ok=True)

#         def makeArrayDatasetFGW(stationFiles):
#             saveFolder = Path(saveFolder)
#             try:
#                 stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx"))

#                 (saveFolder / "temp").mkdir(exist_ok=True)

#                 totalLength = len(stationFiles)
#                 count = 0

#                 for file in tqdm(stationFiles):
#                     count = count + 1
#                     dataFile = preprocess.read_excel_file(file)

#                     if dataFile is None:
#                         continue
#                     logging.info(f"Loaded Excel file: {file}")
#                     wellID = dataFile["ID"].iloc[0]
#                     wellFilePath = file.parents[1] / 'wells.xlsx'
#                     welldataFile = preprocess.read_well_data(wellFilePath)

#                     if welldataFile is None:
#                         continue

#                     locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]

#                     if locationData.empty:
#                         warning_message = f'Well {wellID} is not found'
#                         logging.warning(warning_message)
#                         continue

#                     data = dataFile["Value"]
#                     data = data.values.reshape(data.shape + (1,))

#                     timeVals = dataFile["Date and Time"]
#                     wellVals = np.array([wellID]).astype("str")

#                     dataSet = xr.Dataset(
#                         data_vars={
#                             "Value": ((["time", "well"]), data),
#                         },
#                         coords={
#                             "time": (["time"], timeVals),
#                             "well": (["well"], wellVals),
#                             "lat": (["well"], np.array([locationData["Latitude"].iloc[0]]).astype("float32")),
#                             "lon": (["well"], np.array([locationData["Longitude"].iloc[0]]).astype("float32"))
#                         },
#                     )

#                     if count == 1:
#                         full_ds = dataSet
#                         full_ds['time'] = pd.to_datetime(full_ds['time'])
#                     else:
#                         full_ds = xr.merge([dataSet, full_ds], join="outer", compat="no_conflicts")
#                         full_ds = full_ds.sortby("lat", ascending=False)
#                         full_ds = full_ds.chunk({"well": 1})

#                     full_ds.to_zarr(saveFolder / "gwDepth_array.zarr", mode="w", consolidated=True)
#                     logging.info(f"zarr saved successfully: {file}")

#                 # After processing all files, call makePointsDataset
#                 locationData = makePointsDataset(saveFolder)
#                 return locationData

#             except Exception as e:
#                 error_message = f"Error in preprocessing: {str(e)}"
#                 logging.error(error_message)
#                 print(error_message)
#                 return None

#         def makePointsDataset(station_Files):
#             count = 0
#             obs_wells_data = []

#             # Convert the tuple to a DataFrame
#             locationData = preprocess.locationData
#             locationData_df = pd.DataFrame(locationData, columns=["ID", "Latitude", "Longitude"])
#             print(locationData_df)

#             for index, row in tqdm(locationData_df.iterrows(), total=len(locationData_df), desc="Creating Shapefile"):
#                 stationID = row["ID"]
#                 lat = row["Latitude"]
#                 lon = row["Longitude"]

#                 # Create a Point geometry object
#                 geometry = Point(lon, lat)

#                 # Create a GeoDataFrame for this well
#                 gdf = gpd.GeoDataFrame(
#                     {
#                         "stationID": [stationID],
#                     },
#                     geometry=[geometry],
#                     crs="EPSG:4326"
#                 )

#                 obs_wells_data.append(gdf)
#                 count += 1

#             # Concatenate all the GeoDataFrames into one
#             gdf = pd.concat(obs_wells_data, ignore_index=True)
#             gdf = gpd.GeoDataFrame(gdf, geometry=gdf.geometry, crs="EPSG:4326")

#             # Save the GeoDataFrame to a shapefile
#             (saveFolder / "well_points").mkdir(exist_ok=True)
#             gdf.to_file(saveFolder / "well_points/well_points.shp")
#             return locationData




























#     # @staticmethod
#     # def groundWaterDepthData(groundWaterDepthDataDirectory, saveFolder):
#     #     saveFolder = Path(saveFolder)
#     #     groundWaterDepthDataDirectory = Path(groundWaterDepthDataDirectory)

#     #     saveFolder.mkdir(exist_ok=True)

#     #     try:
#     #         stationFiles = sorted(groundWaterDepthDataDirectory.glob("**/*.xlsx"))
#     #         stationFiles = stationFiles
            

#     #         (saveFolder / "temp").mkdir(exist_ok=True)
            
#     #         totalLength = len(stationFiles)
#     #         count = 0

#     #         for file in tqdm(stationFiles):
#     #             count = count + 1
#     #             dataFile = preprocess.read_excel_file(file)

#     #             if dataFile is None:
#     #                 continue
#     #             logging.info(f"Loaded Excel file: {file}")
#     #             wellID = dataFile["ID"].iloc[0]
#     #             wellFilePath = file.parents[1] / 'wells.xlsx'
#     #             welldataFile = preprocess.read_well_data(wellFilePath)

#     #             if welldataFile is None:
#     #                 continue

#     #             locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]

#     #             if locationData.empty:
#     #                 warning_message = f'Well {wellID} is not found'
#     #                 logging.warning(warning_message)
#     #                 continue
 
#     #             data = dataFile["Value"]
#     #             data = data.values.reshape(data.shape + (1,))
               

#     #             timeVals = dataFile["Date and Time"]
#     #             wellVals = np.array([wellID]).astype("str")
               
                

#     #             dataSet = xr.Dataset(
#     #                 data_vars={
#     #                     "Value": ((["time", "well"]), data),
#     #                 },
#     #                 coords={
#     #                     "time": (["time"], timeVals),
#     #                     "well": (["well"],wellVals),
#     #                     "lat": (["well"], np.array([locationData["Latitude"].iloc[0]]).astype("float32")),
#     #                     "lon": (["well"], np.array([locationData["Longitude"].iloc[0]]).astype("float32"))
#     #                 },
#     #             )

#     #             if count == 1:
#     #                 full_ds = dataSet
#     #                 full_ds['time'] = pd.to_datetime(full_ds['time'])
#     #             else:
#     #                 full_ds = xr.merge([dataSet, full_ds], join="outer", compat="no_conflicts")
#     #                 full_ds = full_ds.sortby("lat", ascending=False)
#     #                 full_ds = full_ds.chunk({"well": 1})
#     #                 print(full_ds)
#     #                 return locationData
                    
                    
                    
#     #             full_ds.to_zarr(saveFolder / "gwDepth_array.zarr", mode="w", consolidated=True)
#     #             logging.info(f"zaar saved sucessfully: {file}")
        

#     #     except Exception as e:
#     #         error_message = f"Error in preprocessing: {str(e)}"
#     #         logging.error(error_message)
#     #         print(error_message)
#     #         return None
#     #             # After processing all files, call makePointsDataset
#     #         locationData = makePointsDataset(saveFolder)
            
#     # def makePointsDataset(locationData, saveFolder):
#     #     saveFolder = Path(saveFolder)
#     #     count = 0
#     #     obs_wells_data = []

#     #     # Convert the tuple to a DataFrame
#     #     locationData = preprocess.locationData
#     #     locationData_df = pd.DataFrame(locationData, columns=["ID", "Latitude", "Longitude"])
#     #     print(locationData_df)

#     #     for index, row in tqdm(locationData_df.iterrows(), total=len(locationData_df), desc="Creating Shapefile"):
#     #         stationID = row["ID"]
#     #         lat = row["Latitude"]
#     #         lon = row["Longitude"]

#     #         # Create a Point geometry object
#     #         geometry = Point(lon, lat)

#     #         # Create a GeoDataFrame for this well
#     #         gdf = gpd.GeoDataFrame(
#     #             {
#     #                 "stationID": [stationID],  
#     #             },
#     #             geometry=[geometry],
#     #             crs="EPSG:4326"
#     #         )

#     #         obs_wells_data.append(gdf)
#     #         count += 1  
#     #     # Concatenate all the GeoDataFrames into one
#     #     gdf = pd.concat(obs_wells_data, ignore_index=True)
#     #     gdf = gpd.GeoDataFrame(gdf, geometry=gdf.geometry, crs="EPSG:4326")  
#     #     # Save the GeoDataFrame to a shapefile
#     #     (saveFolder / "well_points").mkdir(exist_ok=True)
#     #     gdf.to_file(saveFolder / "well_points/well_points.shp")


        


# # Example usage:
# # locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]  # Replace with your data
# # saveFolder = "output_folder_path"
# # makePointsDataset(locationData, saveFolder)

# # Example usage:
# # locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]  # Replace with your data
# # saveFolder = "output_folder_path"
# # makePointsDataset(locationData, saveFolder)


# # Example usage:
# # locationData = welldataFile[welldataFile['ID'].str.match(str(wellID))]  # Replace with your data
# # saveFolder = "output_folder_path"
# # makePointsDataset(locationData, saveFolder)



        
#     # def makePointsDataset(full_ds,saveFolder):
#     #     saveFolder = Path(saveFolder)
#     #     count = 0
#     #     for well in full_ds['well']:
#     #         count = count + 1
#     #         stationInfo = full_ds
#     #         print(stationInfo)
#     #         stationGrid = np.array([stationInfo["well"]])
#     #         stationgGrid = stationgGrid.reshape((1,) + stationgGrid.shape)

#     #         dataSet = pd.DataFrame(
#     #             {
#     #                 "stationID": np.array([stationInfo["well"]]),
#     #                 "lat": np.array([stationInfo["Latitude"]]),
#     #                 "lon": np.array([stationInfo["Longitude"]]),
#     #             }
#     #         )

#     #         if count == 1:
#     #             full_ds = dataSet

#     #         else:
#     #             full_ds = pd.concat([dataSet, full_ds], axis=0, ignore_index=False)

#     #     (saveFolder / "well_points").mkdir(exist_ok=True)
#     #     gdf = gpd.GeoDataFrame(
#     #         full_ds, geometry=gpd.points_from_xy(full_ds.lon, full_ds.lat), crs="EPSG:4326"
#     #     )
#     #     gdf.to_file(saveFolder / "well_points/well_points.shp")
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
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
