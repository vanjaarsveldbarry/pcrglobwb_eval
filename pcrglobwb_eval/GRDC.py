from tqdm import tqdm
from pathlib import Path
import pandas as pd
import warnings
import numpy as np
import xarray as xr
import geopandas as gpd
import pyflwdir
import concurrent.futures
import time
import hydroeval as he
import xarray as xr
from shapely.geometry import Point
import pendulum as pend


def timing_decorator(func):
    '''@timing_decorator ontop of function you want to time'''
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time_min = execution_time / 60
        print(f"\n{func.__name__} took: {execution_time_min:.1f}min / {execution_time:.1f}s to run ")
        return result    
    return wrapper

class GRDC:
    """Preprocess Validation Data for further analysis
    """
    def saveTimeSeries(grdcRawDataDirectory, saveFolder, simDirectory, lddFile):
        """Convert grdc data downloaded from data portal into workable format           
        """   
        saveFolder = Path(saveFolder)
        simDirectory = Path(simDirectory)        
        grdcRawDataDirectory = Path(grdcRawDataDirectory)

        grdcDataFiles = sorted(grdcRawDataDirectory.glob("*.txt"))
        dischargeFiles = sorted(simDirectory.glob(f"**/*netcdf/discharge_dailyTot_output.nc"))
        
        def create_search_criteria(dischargeFiles, no_years=10):
            with xr.open_mfdataset(dischargeFiles,concat_dim='time', engine="netcdf4", 
                                   combine='nested') as ds:
                time_min, time_max = ds.time.min().values, ds.time.max().values
                lat_min, lat_max = ds.lat.min().values.item(), ds.lat.max().values.item()
                lon_min, lon_max = ds.lon.min().values.item(), ds.lon.max().values.item()
                return {'time_min': time_min, 'time_max': time_max, 
                        'lat_min': lat_min,'lat_max': lat_max,
                        'lon_min': lon_min,'lon_max': lon_max,
                        'no_years': no_years}

        def checkOverlap(search_criteria, result, catchmentSize = 5.0, yearsOverlap=3):
            # Check if the time ranges overlap
            if result['ts_end'] < search_criteria['time_min'] or result['ts_start'] > search_criteria['time_max']:
                return False
            
            # Check if the spatial ranges overlap
            if search_criteria['lat_min'] < result['latitude'] < search_criteria['lat_max']:
                if search_criteria['lon_min'] < result['longitude'] < search_criteria['lon_max']:
                    if result['cat_area'] >= catchmentSize:
                        ts_end = pend.from_format(result['ts_end'].strftime("%Y-%m-%d"), 'YYYY-MM-DD')
                        ts_start = pend.from_format(result['ts_start'].strftime("%Y-%m-%d"), 'YYYY-MM-DD')
                        time_min = pend.from_format(search_criteria['time_min'].astype(str)[:10], 'YYYY-MM-DD')
                        time_max = pend.from_format(search_criteria['time_max'].astype(str)[:10], 'YYYY-MM-DD')
                        
                        if ts_end < time_min or ts_start > time_max:
                            return False
                        # Check if the difference between the start and end dates is at least 3 years
                        if ts_start > time_min:
                            latest_start = ts_start
                        else:
                            latest_start = time_min
                        if ts_end < time_max:
                            earliest_end = ts_end
                        else:
                            earliest_end = time_max

                        overlap = (earliest_end - latest_start).days if latest_start <= earliest_end else 0
                        if  overlap >= yearsOverlap*365:
                            return True
                        else:
                            return False
                    else:
                        return False
                else:
                    return False    
            else:
                return False     
            
        def search_grdc_stations(f, search_criteria, encoding='ISO-8859-1'):
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

                if 'Latitude' in split_line:
                    lat_grdc = float(line.split(":")[-1].strip())
                    props['latitude'] = float(lat_grdc)
                
                if 'Longitude' in split_line:
                    lon_grdc = line.split(":")[-1].strip()
                    props['longitude'] = float(lon_grdc)
                    
                if 'GRDC-No.' in split_line:
                    grdc_no = line.split(":")[-1].strip()
                    props['grdc_no'] = int(grdc_no)

                if 'Catchment area' in split_line:
                    cat_area = line.split(":")[-1].strip()
                    if cat_area == '':
                        cat_area = 0.0
                    props['cat_area'] = float(cat_area)

                if 'No. of years' in split_line:
                    no_years = line.split(":")[-1].strip()
                    props['no_years'] = int(no_years)
                
                if 'Time series' in split_line:
                    ts_time = line.split(":")[-1].strip()
                    ts_start = ts_time.split(" - ")[0].strip()
                    ts_end = ts_time.split(" - ")[1].strip()
                    props['ts_start'] = pd.to_datetime(ts_start)
                    props['ts_end'] = pd.to_datetime(ts_end) 

                # break loop to save time and not read all observed values   
                if i > 25: break
                    
            # close file        
            f.close()
            
            if checkOverlap(search_criteria, props) == True: 
                return props['grdc_no']
            else: 
                return None

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

                if 'Latitude' in split_line:
                    lat_grdc = float(line.split(":")[-1].strip())
                    props['latitude'] = float(lat_grdc)
                
                if 'Longitude' in split_line:
                    lon_grdc = line.split(":")[-1].strip()
                    props['longitude'] = float(lon_grdc)
                    
                if 'GRDC-No.' in split_line:
                    grdc_no = line.split(":")[-1].strip()
                    props['grdc_no'] = int(grdc_no)

                # if 'Station' in split_line:
                #     station_grdc = line.split(":")[-1].strip()
                #     props['station'] = str(station_grdc)


                if 'Catchment area' in split_line:
                    cat_area = line.split(":")[-1].strip()
                    if cat_area == '':
                        cat_area = 0.0
                    props['cat_area'] = float(cat_area)

                if 'No. of years' in split_line:
                    no_years = line.split(":")[-1].strip()
                    props['no_years'] = int(no_years)
                
                if 'Time series' in split_line:
                    ts_time = line.split(":")[-1].strip()
                    ts_start = ts_time.split(" - ")[0].strip()
                    ts_end = ts_time.split(" - ")[1].strip()
                    props['ts_start'] = pd.to_datetime(ts_start)
                    props['ts_end'] = pd.to_datetime(ts_end) 

                # break loop to save time and not read all observed values   
                if i > 25: break
                    
            # close file        
            f.close()
            return props

        def get_grdc_station_values(f, var_name=None, col_name=' Value', remove_mv=True, mv_val=-999, encoding='ISO-8859-1', verbose=False):
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
                                data=df[var_name]).rename(columns={var_name: 'grdc_discharge'})
            if remove_mv == True:
                df_out.replace(mv_val, np.nan, inplace=True)

            return df_out

        def correctLatLons(grdcStations, lddFile, dischargeFiles, search_criteria):
            def adjustPoint(props, upstream_area, res):
                window = 5 * res
                lonMin, lonMax = (props['longitude'] - window), (props['longitude'] + window)
                latMin , latMax = (props['latitude'] - window), (props['latitude'] + window)
                sub_upstream_area = upstream_area.sel(lon=slice(lonMin, lonMax), lat=slice(latMax, latMin)).compute()
                diff = np.absolute(sub_upstream_area.upstream_area - props['cat_area'])
                xyID = diff.argmin(dim=["lat", "lon"])
                sub_upstream_area = sub_upstream_area.upstream_area.isel(lat=np.array([xyID["lat"].values]), 
                                                                lon=np.array([xyID["lon"].values])).squeeze(["lat", "lon"])
                newLon = sub_upstream_area.lon.values.item()
                newLat = sub_upstream_area.lat.values.item()
                
                # Update props with new latitude and longitude values
                props['latitude'] = newLat
                props['longitude'] = newLon
                return props        
            def quick_pass_check(props ,sim_ds):
                ds = sim_ds.sel(lon=props['longitude'], lat=props['latitude'], method="nearest").values
                if np.isnan(ds):
                    return False
                else:
                    return True
            sim_ds = xr.open_dataset(dischargeFiles[0], decode_coords=False, decode_times=False, chunks='auto', engine='h5netcdf')
            
            sim_ds = sim_ds.isel(time=0).discharge
            res = sim_ds['lat'].values[0] - sim_ds['lat'].values[1]
            if res < 0.5:  
                with xr.open_dataset(lddFile, chunks='auto') as upstream_area:
                    upstream_area = upstream_area.sel(lon=slice(search_criteria['lon_min'], search_criteria['lon_max']), 
                                        lat=slice(search_criteria['lat_max'], search_criteria['lat_min']))
                    upstream_area= upstream_area.rename_vars({var: "upstream_area" for var in upstream_area.data_vars})
                    template_array = upstream_area.copy(deep=False)
                    upstream_area = upstream_area.rename({'lat': 'latitude', 'lon':'longitude'})
                    upstream_area = upstream_area.rio.write_crs("epsg:4326")
                    transform = upstream_area.rio.transform()
                    upstream_area = upstream_area.upstream_area.values.astype(np.uint8)
                    upstream_area[upstream_area == 0] = 255
                    upstream_area = pyflwdir.from_array(upstream_area, ftype='ldd', transform=transform, latlon=True, cache=True)
                    upstream_area = upstream_area.upstream_area(unit='km2')
                    upstream_area = template_array.copy(data={"upstream_area": upstream_area})
                filter_grdcStations = []
                for station in tqdm(grdcStations, desc="Correcting Coordinates", disable=False):
                    matching_files = [str(f) for f in grdcDataFiles if str(station) in str(f)][0]
                    props = get_grdc_station_properties(matching_files)
                    if quick_pass_check(props, sim_ds) == True:
                        props = adjustPoint(props, upstream_area, res)
                        filter_grdcStations.append(props)
                sim_ds.close()        
                return filter_grdcStations    

            else:
                filter_grdcStations = []
                for station in tqdm(grdcStations, desc="Correcting Coordinates", disable=False):
                    matching_files = [str(f) for f in grdcDataFiles if str(station) in str(f)][0]
                    props = get_grdc_station_properties(matching_files)
                    filter_grdcStations.append(props)
                sim_ds.close()        
                return filter_grdcStations    
        
        def getIndices(dischargeFiles, lats, lons):
            with xr.open_mfdataset(dischargeFiles, combine="nested",concat_dim=["time"], decode_coords=False, 
                            combine_attrs="drop", decode_times=False, engine='h5netcdf',  decode_cf=False, parallel=True, chunks={'time':-1, 'lat':10, 'lon':10}) as sim_ds:
                lon_indices = []
                for lon in lons:
                    lon_diff = np.abs(sim_ds.coords["lon"].values - lon)
                    lon_index = np.where(lon_diff == np.min(lon_diff))[0][0]
                    lon_indices.append(lon_index)
                lat_indices = []
                for lat in lats:
                    lat_diff = np.abs(sim_ds.coords["lat"].values - lat)
                    lat_index = np.where(lat_diff == np.min(lat_diff))[0][0]
                    lat_indices.append(lat_index)
                if lat_indices == [] or lon_indices == []:
                    return None, None
                lat_indices = xr.DataArray(lat_indices, dims='station', coords={"station": [station['grdc_no'] for station in grdcStations]}).astype('int32')
                lon_indices = xr.DataArray(lon_indices, dims='station', coords={"station": [station['grdc_no'] for station in grdcStations]}).astype('int32')
                return lat_indices, lon_indices

        def read_netcdfs(dischargeFiles, lon_indices, lat_indices, step_size=1):
            data = xr.open_mfdataset(dischargeFiles, decode_cf=False,
                                   concat_dim='time', engine="h5netcdf", 
                                   combine='nested', combine_attrs='drop',
                                   join='override', parallel=True,
                                   mask_and_scale=False, 
                                   chunks=None)['discharge']
            timeSteps = [i for i in range(0, data['time'].size+1, step_size)]
            for timeStep in tqdm(timeSteps, desc="Reading Points", disable=False):
                data_sub = data.isel(time=slice(timeStep,timeStep+step_size)).load().isel(lon=lon_indices, lat=lat_indices)
                if timeStep == 0: dsList = [data_sub]
                else: dsList.append(data_sub)
            selectedData = xr.concat(dsList, 'time')
            return selectedData
        
        def create_gpd(ds):
            df = pd.DataFrame({'lat': ds.lat.values, 'lon': ds.lon.values, 'station': ds.station.values})
            geometry = [Point(xy) for xy in zip(df.lon, df.lat)]
            gdf = gpd.GeoDataFrame(df, geometry=geometry)
            gdf = gdf.drop(['lat', 'lon'], axis=1)
            return gdf.to_file(saveFolder / f"points/points.shp")

        search_criteria = create_search_criteria(dischargeFiles)
        grdcStations = []
        with concurrent.futures.ThreadPoolExecutor(25) as executor:
            futures = [executor.submit(search_grdc_stations, file, search_criteria) for file in grdcDataFiles]
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Filter grdc stations", disable=False):
                result = future.result()
                if result is not None:
                    grdcStations.append(result)
        
        grdcStations = correctLatLons(grdcStations, lddFile, dischargeFiles, search_criteria)
        lat_indices, lon_indices = getIndices(dischargeFiles,grdcStations)

        if lat_indices is not None and lon_indices is not None:
            (saveFolder / f"timeseries").mkdir(exist_ok=True, parents=True)
            (saveFolder / f"scores").mkdir(exist_ok=True, parents=True) 
            (saveFolder / f"points").mkdir(exist_ok=True, parents=True) 
            
            sim_ds = read_netcdfs(dischargeFiles, lon_indices, lat_indices)
            create_gpd(sim_ds)
            for station in tqdm(grdcStations, desc="Processing grdc stations", disable=False):
                matching_files = [str(f) for f in grdcDataFiles if str(station['grdc_no']) in str(f)][0]
                df_values = get_grdc_station_values(matching_files)
                df_values = df_values.loc[(df_values.index >= search_criteria['time_min']) & (df_values.index <= search_criteria['time_max'])]
                time_range = pd.date_range(start=search_criteria['time_min'], end=search_criteria['time_max'], freq='1D')
                df_values = df_values.reindex(time_range)
                df_values['sim_discharge']=sim_ds.sel(station=station['grdc_no']).values
                first_non_missing, last_non_missing = df_values['grdc_discharge'].first_valid_index(), df_values['grdc_discharge'].last_valid_index()
                df_values = df_values.loc[first_non_missing:last_non_missing]
                non_nan_count = df_values['grdc_discharge'].count()
                if non_nan_count > 1095:
                    df_values.to_parquet(saveFolder / f"timeseries/{station['grdc_no']}.parquet")
            sim_ds.close()

    def scoreTimeSeries(saveFolder):
        saveFolder = Path(saveFolder)
        timeSeriesFiles = sorted(saveFolder.glob(f"timeseries/*.parquet"))
        if len( timeSeriesFiles) != 0:
            df = pd.DataFrame(columns=['station', 'KGE', 'KGE_r', 'KGE_alpha', 'KGE_beta', 'KGE_bnd'])
            for file in tqdm(timeSeriesFiles, desc="Scoring Stations Daily", disable=False):
                station=int(file.stem)
                file = pd.read_parquet(file)
                kge, r, alpha, beta = he.evaluator(he.kge, file["sim_discharge"], file["grdc_discharge"])
                kge_bounded = he.evaluator(he.kge_c2m, file["sim_discharge"], file["grdc_discharge"])

                if df.empty:
                    df = pd.DataFrame({'station': [station], 
                                       'KGE': kge,
                                       'KGE_r': r, 
                                       'KGE_alpha': alpha, 
                                       'KGE_beta': beta,
                                       'KGE_bnd': kge_bounded})
                else:
                    df = pd.concat([df, pd.DataFrame({'station': [station], 
                                                      'KGE': kge,
                                                      'KGE_r': r, 
                                                      'KGE_alpha': alpha,
                                                      'KGE_beta': beta,
                                                      'KGE_bnd': kge_bounded})], ignore_index=True)
            scores = gpd.read_file(saveFolder / "points/points.shp")
            scores = scores.merge(df, on="station")
            scores.to_file(saveFolder / f"scores/day_scores.shp")
            
            df = pd.DataFrame(columns=['station', 'KGE', 'KGE_r', 'KGE_alpha', 'KGE_beta', 'KGE_bnd'])
            for file in tqdm(timeSeriesFiles, desc="Scoring Stations Monthly", disable=False):
                station=int(file.stem)
                file = pd.read_parquet(file)
                file = file.resample('M').mean()
                kge, r, alpha, beta = he.evaluator(he.kge, file["sim_discharge"], file["grdc_discharge"])
                if df.empty:
                    df = pd.DataFrame({'station': [station], 
                                       'KGE': kge,
                                       'KGE_r': r, 
                                       'KGE_alpha': alpha, 
                                       'KGE_beta': beta,
                                       'KGE_bnd': kge_bounded})
                else:
                    df = pd.concat([df, pd.DataFrame({'station': [station], 
                                                      'KGE': kge,
                                                      'KGE_r': r, 
                                                      'KGE_alpha': alpha,
                                                      'KGE_beta': beta,
                                                      'KGE_bnd': kge_bounded})], ignore_index=True)
            scores = gpd.read_file(saveFolder / "points/points.shp")
            scores = scores.merge(df, on="station")
            scores.to_file(saveFolder / "scores/month_scores.shp")