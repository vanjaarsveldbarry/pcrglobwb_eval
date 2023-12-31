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
import concurrent.futures
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
    
    
    return wrapper
   

class preprocess:
    """Preprocess Validation Data for further analysis
    """
    def grdcData(grdcRawDataDirectory, saveFolder, lddFile):
        """Convert grdc data downloaded from data portal to:
            - GRDC_array.zarr
            - GRDC_points.shp
        """        
        
        @timing_decorator   
        def makeArrayDataset(grdcDataFiless):
            
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
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=80) as executor:
                futures = [executor.submit(get_grdc_station_properties, file) for file in grdcDataFiles]
                dataSetList = []
                for future in futures: dataSetList.append(future.result())
             
            startVals = [d.get('ts_start') for d in dataSetList]
            endVals = [d.get('ts_end') for d in dataSetList]
            timeVals = pd.date_range(min(startVals), max(endVals), freq='D')
            stationNumber = [d.get('grdc_no') for d in dataSetList]
            latVals = [d.get('latitude') for d in dataSetList]
            lonVals = [d.get('longitude') for d in dataSetList]
            cat_area = [d.get('cat_area') for d in dataSetList]
            #TODO reincorporate the name
            # stationName = [d.get('station') for d in dataSetList]
            
            dataSet = xr.Dataset(
                        data_vars={
                            "grdc_discharge": ((["station", "time"]), np.empty((len(stationNumber), timeVals.shape[0]))) 
                        },
                        coords={
                            "time": (["time"], timeVals),
                            "station": (["station"], np.array(stationNumber).astype("int32")),
                            "lat": (["station"], np.array(latVals).astype("float32")),
                            "lon": (["station"], np.array(lonVals).astype("float32")),
                            "ts_start": (["station"], np.array(startVals)),	
                            "ts_end": (["station"], np.array(endVals)),
                            "cat_area": (["station"], np.array(cat_area).astype("float32")),
                            #TODO include name
                            # "name": (["station"], np.array([f"{stationInfo['station']} (No. {stationInfo['grdc_no']})"]).astype("str")),
                        })
            dataSet.to_zarr(saveFolder / "GRDC_array.zarr", mode="w", compute=False,consolidated=True) 
            
            def process_grdc_file(file):
                
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

                    return df 

                stationInfo = get_grdc_station_properties(file)            
                data = get_grdc_station_values(file, var_name="OBS", col_name=" Value")

                dataSet = xr.Dataset(
                    data_vars={"grdc_discharge": ((["station", "time"]), data["OBS"].values[np.newaxis, ...].astype("float32"))},
                    coords={
                        "time": (["time"], data["OBS"].index.values),
                        "station": (["station"], np.array([stationInfo["grdc_no"]]).astype("int32")),
                        "lat": (["station"], np.array([stationInfo["latitude"]]).astype("float32")),
                        "lon": (["station"], np.array([stationInfo["longitude"]]).astype("float32")),
                        "ts_start": (["station"], np.array([stationInfo["ts_start"]])),
                        "ts_end": (["station"], np.array([stationInfo["ts_end"]])),
                        "cat_area": (["station"], np.array([stationInfo["cat_area"]]).astype("float32")),
                        #TODO INCLUDE NAME
                        # "name": (["station"], np.array([f"{stationInfo['station']} (No. {stationInfo['grdc_no']})"]).astype("str")),
                    })
                dataSet = dataSet.reindex(time=pd.date_range(stationInfo["ts_start"], stationInfo["ts_end"], freq='D'))
                return dataSet

            def write_chunk_file(file):
                data = process_grdc_file(file)
                stationIndex = np.where(dataSet.station.values == data.station)[0][0]
                timeStartIndex = np.where(dataSet.time.values == data['ts_start'])[0][0]
                timeEndIndex = np.where(dataSet.time.values == data['ts_end'])[0][0]
                data.to_zarr(saveFolder / "GRDC_array.zarr", region={"time": slice(timeStartIndex, timeEndIndex+1), "station": slice(stationIndex,stationIndex+1)})

            with concurrent.futures.ProcessPoolExecutor(max_workers=80) as executor:
                futures = []
                for file in grdcDataFiles: futures.append(executor.submit(write_chunk_file, file))   
                for future in tqdm(concurrent.futures.as_completed(futures)): pass        

        @timing_decorator   
        def makePointsDataset(saveFolder):
                # read array then pullout ID, lat, lon.
                (saveFolder / "GRDC_points").mkdir(exist_ok=True)
                ds = xr.open_zarr(saveFolder / "GRDC_array.zarr")
                df = pd.DataFrame({'lat': ds.lat.values, 'lon': ds.lon.values, 'station': ds.station.values})
                gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4326")
                gdf.to_file(saveFolder / "GRDC_points/GRDC_points.shp")

        @timing_decorator
        def makeCatchmentAreaMap(lddFile):
            """Make a catchment area map based on the ldd file
            """
            ldd_ds = xr.open_dataset(lddFile, chunks='auto')
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
        
        grdcDataFiles = sorted(grdcRawDataDirectory.glob("*.txt"))
        makeArrayDataset(grdcDataFiles)
        makePointsDataset(saveFolder)
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