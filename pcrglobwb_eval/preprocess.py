from tqdm import tqdm
from pathlib import Path
import pandas as pd
import warnings
import numpy as np
import xarray as xr
import shutil
import geopandas as gpd
from functools import partial
import pyflwdir

import time
def timing_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"{func.__name__} took {execution_time:.6f} seconds to run")
        return result
    return wrapper
# @timing_decorator

class preprocess:
    """Preprocess Validation Data for further analysis
    """
    @timing_decorator
    def grdcData(grdcRawDataDirectory, saveFolder, lddFile):
        """Convert grdc data downloaded from data portal to:
            - GRDC_array.zarr
            - GRDC_points.shp
        """

        def get_grdc_station_properties(f, encoding='ISO-8859-1'):
            """Retrieves GRDC station properties from txt-file. Creates and returns header from those properties as well as a dictionary containt station name, lat, and lon info.

            Args:
                encoding (str, optional): encoding of GDRC files. Defaults to 'ISO-8859-1'.

            Returns:
                dict: dictionary containing properties.
            """

            # open file
            f = open(f, encoding=encoding)

            props = dict()
            
            # go through lines in file
            for i, line in enumerate(f):

                split_line = line.split(":")[0]

                if 'GRDC-No.' in split_line:
                    grdc_no = line.split(":")[-1].strip()
                    props['grdc_no'] = int(grdc_no)

                if 'Station' in split_line:
                    station_grdc = line.split(":")[-1].strip()
                    props['station'] = str(station_grdc)

                if 'Latitude' in split_line:
                    lat_grdc = line.split(":")[-1].strip()
                    props['latitude'] = float(lat_grdc)

                if 'Longitude' in split_line:
                    lon_grdc = line.split(":")[-1].strip()
                    props['longitude'] = float(lon_grdc)

                if 'Catchment area' in split_line:
                    cat_area = line.split(":")[-1].strip()
                    if cat_area == '':
                        cat_area = 0.0
                    props['cat_area'] = float(cat_area)

                if 'Time series' in split_line:
                    ts_time = line.split(":")[-1].strip()
                    ts_start = ts_time.split(" - ")[0].strip()
                    ts_end = ts_time.split(" - ")[1].strip()
                    props['ts_start'] = pd.to_datetime(ts_start)
                    props['ts_end'] = pd.to_datetime(ts_end)

                if 'No. of years' in split_line:
                    no_years = line.split(":")[-1].strip()
                    props['no_years'] = int(no_years)

                # break loop to save time and not read all observed values   
                if i > 25:
                    break
                    
            # close file        
            f.close()

            if 'grdc_no' not in props.keys(): warnings.warn('WARNING -- no "GRDC-No." information found in file.')
            if 'station' not in props.keys(): warnings.warn('WARNING -- no "Station" information found in file.')
            if 'latitude' not in props.keys(): warnings.warn('WARNING -- no "Latitude" information found in file.')
            if 'longitude' not in props.keys(): warnings.warn('WARNING -- no "Longitude" information found in file.')
            if 'cat_area' not in props.keys(): warnings.warn('WARNING -- no "Catchment area" information found in file.')
            if 'ts_start' not in props.keys(): warnings.warn('WARNING -- no start date of timeseries found in file.')
            if 'ts_end' not in props.keys(): warnings.warn('WARNING -- no end date of timeseries found in file.')
            if 'no_years' not in props.keys(): warnings.warn('WARNING -- no "No. of years" information found in file.')

            return props

        def get_grdc_station_values(f, var_name=None, col_name=' Value', remove_mv=True, mv_val=-999, encoding='ISO-8859-1', verbose=False) -> pd.DataFrame:
            """Reads (discharge-)values of GRDC station from txt-file and returns them as dataframe. 
            Creates a pandas dataframe with a user-specified column header for values instead of default ' Values' header name. 
            Possible to remove possible missing values in the timeseries and plot the resulting series.

            Args:
                var_name (str): user-specified variable name to be given to column. If None, col_name is used. Default to 'None'.
                col_name (str, optional): name of column in GRDC-file to be read. Defaults to ' Value'.
                remove_mv (bool, optional): whether or not remove missing values in timeseries. Defaults to True.
                mv_val (int, optional): missing value in timeseries. Defaults to -999.
                encoding (str, optional): encoding of GDRC files. Defaults to 'ISO-8859-1'.
                verbose (bool, optional): whether or not to show more info. Defaults to False.
            
            Returns:
                pd.DataFrame: dataframe containing observational data.
            """

            # open file
            f = open(f, encoding=encoding)

            # find line in file in which meta-data starts and observational record starts
            # for i, line in enumerate(f):
            #     if '#' in line:
            #         pass
            #     else:
            #         stopline = i-1
            #         break

            df = pd.read_csv(f, header=0, comment='#', sep=';', encoding=encoding)

            # if var_name is specified, use it
            if var_name != None:
                var_name = var_name
            else:
                var_name = str(col_name)

            df[var_name] = df[str(col_name)].copy()

            try: 
                if verbose: click.echo('VERBOSE -- reading column {}'.format(col_name))
                df[var_name] = df[str(col_name)].copy()

            except:
                if col_name == ' Value':
                    raise ValueError(f'ERROR: column "{col_name}" - which is also the fall back option - cannot be found in file {f}')
                else:
                    warnings.warn('WARNING: column {} not found, falling back to column Value'.format(col_name))
                    df[var_name] = df[' Value'].copy()
                    del df[' Value']

            #To deal with instances where YYYY-MM-DD is YYY-MM-00, pd.to_datetime() throws error for 00.
            if int(df['YYYY-MM-DD'][0][-1:]) == 0: 
                df['YYYY-MM-DD'] =  df['YYYY-MM-DD'].str.replace('-00', '-01')
                
            df['date'] = pd.to_datetime(df['YYYY-MM-DD'])
            df.set_index(df['date'], inplace=True)

            df_out = pd.DataFrame(index=df.index,
                                data=df[var_name])
            if remove_mv == True:
                df_out.replace(mv_val, np.nan, inplace=True)
            # if (pd.infer_freq(df_out.index) == 'M') or (pd.infer_freq(df_out.index) == 'MS'):
            #     # if verbose: print('changing index strftime to %Y-%m')
            #     df_out.index = df_out.index.strftime('%Y-%m')
            # self.df = df_out

            return df 
        @timing_decorator   
        def makeArrayDataset(dayFiles):
            
            (saveFolder / "temp").mkdir(exist_ok=True)
            totalLength = len(dayFiles)
            count = 0
            for file in tqdm(dayFiles):
                count = count + 1
                stationInfo = get_grdc_station_properties(file)            
                data = get_grdc_station_values(file, var_name="OBS", col_name=" Value")
                data = data["OBS"]

                timeVals = data.index.values
                data = data.values
                data = data.reshape(data.shape + (1,))


                dataSet = xr.Dataset(
                    data_vars={
                        "grdc_discharge": ((["time", "station"]), data),
                    },
                    coords={
                        "time": (["time"], timeVals),
                        "station": (["station"], np.array([stationInfo["grdc_no"]]).astype("int32")),
                        "lat": (["station"], np.array([stationInfo["latitude"]]).astype("float32")),
                        "lon": (["station"], np.array([stationInfo["longitude"]]).astype("float32")),
                        "ts_start": (["station"], np.array([stationInfo["ts_start"]])),
                        "ts_end": (["station"], np.array([stationInfo["ts_end"]])),
                        "cat_area": (["station"], np.array([stationInfo["cat_area"]]).astype("int32")),
                        "name": (
                            ["station"], np.array([f"{stationInfo['station']} (No. {stationInfo['grdc_no']})"]).astype("str"),
                        ),
                    },
                )
                if count == 1:
                    full_ds = dataSet

                else:
                    dataSet, full_ds = xr.align(dataSet, full_ds, exclude=["lat", "lon"], join="outer")

                    full_ds = xr.merge([dataSet, full_ds], join="outer", compat="no_conflicts")
                    full_ds = full_ds.sortby("lat", ascending=False)

                if not count % 10 or count == totalLength:
                    full_ds = full_ds.isel(station=slice(None, -1)).astype("float32").compute()
                    full_ds.to_zarr(saveFolder / f"temp/test{count}.zarr", mode="w")
                    full_ds = full_ds.isel(station=slice(-2, -1))
            
            fullFiles = sorted((saveFolder / f"temp").glob("*/"))
            full_ds = xr.open_mfdataset(fullFiles, chunks=None, parallel=True, engine="zarr", combine="nested", concat_dim="station")
            full_ds = full_ds.drop_duplicates(dim="station")
            full_ds = full_ds.reset_encoding()
            full_ds.to_zarr(saveFolder / "GRDC_array.zarr", mode="w", consolidated=True)
            shutil.rmtree(saveFolder / "temp")
            
        @timing_decorator   
        def makePointsDataset(dayFiles):
                count = 0
                for file in tqdm(dayFiles):
                    count = count + 1
                    stationInfo = get_grdc_station_properties(file)            
                    stationgGrid = np.array([stationInfo["grdc_no"]])
                    stationgGrid = stationgGrid.reshape((1,) + stationgGrid.shape)

                    dataSet = pd.DataFrame(
                        {"stationID": np.array([stationInfo["grdc_no"]]),
                         "lat": np.array([stationInfo["latitude"]]),
                         "lon": np.array([stationInfo["longitude"]]),
                        }
                    )

                    if count == 1:
                        full_ds = dataSet

                    else:
                        full_ds = pd.concat([dataSet, full_ds], axis=0, ignore_index=False)

                (saveFolder / "GRDC_points").mkdir(exist_ok=True)
                gdf = gpd.GeoDataFrame(full_ds, geometry=gpd.points_from_xy(full_ds.lon, full_ds.lat), crs="EPSG:4326")
                gdf.to_file(saveFolder / "GRDC_points/GRDC_points.shp")

        @timing_decorator
        def makeCatchmentAreaMap(lddFile):
            ldd_ds = xr.open_dataset(lddFile)
            ldd_ds = ldd_ds.rename({'lat': 'latitude', 'lon':'longitude'})
            ldd = ldd_ds.where(ldd_ds.Band1!=0., 255.)
            ldd = ldd.rio.write_crs("epsg:4326")
            transform = ldd.rio.transform()
            ldd = ldd.Band1.values.astype(np.uint8)
            upstream_area = pyflwdir.from_array(ldd, ftype='ldd', transform=transform, latlon=True, cache=True)
            upstream_area = upstream_area.upstream_area(unit='km2')
            upstream_area = ldd_ds.rename({'Band1': 'upstream_area'}).copy(data={"upstream_area": upstream_area})
            upstream_area.attrs = {}
            ldd_ds = ldd_ds.rename({'latitude':'lat', 'longitude':'lon'})
            upstream_area.to_zarr(saveFolder / "GRDC_upstreamArea.zarr", mode="w", consolidated=True)
            
        saveFolder = Path(saveFolder)
        saveFolder.mkdir(exist_ok=True)
        
        grdcRawDataDirectory = Path(grdcRawDataDirectory)
        
        dayFiles = sorted(grdcRawDataDirectory.glob("*.txt"))
        makeArrayDataset(dayFiles)
        makePointsDataset(dayFiles)
        makeCatchmentAreaMap(lddFile)
        
    @timing_decorator
    def snowData(snowRawDataDirectory, saveFolder):
        """Convert raw snow data to processed form
            - snow_array.zarr
        """
        
        def makeArrayDatasetSnow(snowFiles):

            def _preprocess(ds):
                ds = ds.sortby('lat', ascending=False)
                ds = ds.where(~ds.isin([255,254,253,252,210,205,206])) #remove error cells
                ds = xr.where(ds > 0, 1., ds.scfg)                #convert snow cover percentage to binary & galcier, snow and ice
                return ds[['scfg']]
            
            partial_func = partial(_preprocess)
        
            ds = xr.open_mfdataset(snowFiles, chunks=({'time': -1, 'lat': 1000, 'lon': 1000}), parallel=True, concat_dim="time", combine="nested", data_vars='minimal', 
                                   coords='minimal', engine='h5netcdf', compat='override', preprocess=partial_func)
            ds = ds.to_zarr(saveFolder / 'snowCover.zarr', mode='w')
        
        saveFolder = Path(saveFolder)
        saveFolder.mkdir(exist_ok=True)
        snowRawDataDirectory = Path(snowRawDataDirectory)
        snowRawDataDirectory = snowRawDataDirectory
        snowFiles = sorted(snowRawDataDirectory.glob("**/*.nc"))
        makeArrayDatasetSnow(snowFiles)
        
    # def graceData(graceRawDataDirectory, saveFolder)