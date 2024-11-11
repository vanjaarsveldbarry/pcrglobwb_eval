from tqdm import tqdm
from pathlib import Path
import numpy as np
import xarray as xr
import time
import xarray as xr
import hydroeval as he
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import pendulum as pend
pd.set_option('display.max_columns', None, 'display.max_rows', None)


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

class soilmoisture:
    def saveTimeSeries(fluxRawDataDirectory, layerDepthFile, saveFolder, simDirectory):
        """Convert grdc data downloaded from data portal into workable format           
        """   
        saveFolder = Path(saveFolder)
        simDirectory = Path(simDirectory)        
        fluxRawDataDirectory = Path(fluxRawDataDirectory)

        fluxMetaFile = fluxRawDataDirectory / 'meta.csv'
        fluxFiles = sorted(fluxRawDataDirectory.glob('data/**/*.stm'))

        soilmositureFiles = sorted(simDirectory.glob(f"**/*netcdf/satDegUpp_dailyTot_output.nc"))
        def create_search_criteria(soilmositureFiles, no_years=2):
            with xr.open_mfdataset(soilmositureFiles, chunks={}, concat_dim='time', engine="netcdf4", combine='nested') as ds:
                # ds = ds.isel(time=slice(0,365*no_years))
                time_min, time_max = ds.time.min().dt.strftime('%Y-%m-%d').item(), ds.time.max().dt.strftime('%Y-%m-%d').item()
                lat_min, lat_max = ds.lat.min().values.item(), ds.lat.max().values.item()
                lon_min, lon_max = ds.lon.min().values.item(), ds.lon.max().values.item()
                return {'time_min': time_min, 'time_max': time_max, 
                        'lat_min': lat_min,'lat_max': lat_max,
                        'lon_min': lon_min,'lon_max': lon_max,
                        'no_years': no_years}
                
        def findFluxSites(fluxMetaFile, search_criteria):
            def checkOverlap(search_criteria, result, yearsOverlap=3):
                #Check if the time ranges overlap
                try:
                    ts_end = pend.from_format(result['ts_end'], 'YYYY-MM-DD')
                    ts_start = pend.from_format(result['ts_start'], 'YYYY-MM-DD')
                    time_min = pend.from_format(search_criteria['time_min'], 'YYYY-MM-DD')
                    time_max = pend.from_format(search_criteria['time_max'], 'YYYY-MM-DD')
                    if ts_end < time_min or ts_start > time_max:
                        return False
                    
                    # Check if the spatial ranges overlap
                    if search_criteria['lat_min'] < result['latitude'] < search_criteria['lat_max']:
                        if search_criteria['lon_min'] < result['longitude'] < search_criteria['lon_max']:
                            # Check if the difference between the start and end dates is at least 3 years
                            if ts_start > time_min:
                                latest_start = ts_start
                            else:
                                latest_start = time_min
                            if ts_end < time_max:
                                earliest_end = ts_end
                            else:
                                earliest_end = time_max

                            # Calculate the overlapping days
                            overlap = (earliest_end - latest_start).days if latest_start <= earliest_end else 0
                            if  overlap >= yearsOverlap*365:
                                return True
                            else:
                                return False
                        else:
                            return False
                    else:
                        return False 
                except:
                    pass
                

            meta_df = pd.read_csv(fluxMetaFile)
            df = pd.DataFrame(columns=['latitude', 'longitude', 'ts_start', 'ts_end'])
            for index, row in meta_df.iterrows():
                result = {'id': row['id'],
                          'name': row['name'],
                          'latitude': row['lat'], 'longitude': row['lon'],
                          'ts_start': row['start-date'], 
                          'depth_end': row['depth_end'],
                          'ts_end': row['end-date'],
                          'file_name': row['file_name']}
                
                if checkOverlap(search_criteria, result) == True:
                    df = pd.concat([df, pd.DataFrame(
                        {'id': [result['id']],
                         'name': [result['name']],
                         'latitude': [result['latitude']], 
                         'longitude':[result['longitude']],
                         'depth_end': result['depth_end'],
                         'ts_start': [result['ts_start']], 
                         'ts_end': [result['ts_end']], 
                         'file_name': row['file_name']})], ignore_index=True)
            return df
        
        def getIndices(soilmositureFiles, metaInfo):
            lats = metaInfo['latitude'].values
            lons = metaInfo['longitude'].values
            with xr.open_mfdataset(soilmositureFiles, combine="nested",concat_dim=["time"], decode_coords=False, 
                            combine_attrs="drop", decode_times=False, engine='h5netcdf',  decode_cf=False, parallel=True, chunks={}) as sim_ds:

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
                lat_indices = xr.DataArray(lat_indices, dims='station', coords={"station": metaInfo['id'].values}).astype('int32')
                lon_indices = xr.DataArray(lon_indices, dims='station', coords={"station": metaInfo['id'].values}).astype('int32')
                ds = sim_ds.isel(time=0).load().isel(lon=lon_indices, lat=lat_indices)
                ds = ds.where(ds < 1e+20)
                ds = ds.dropna(dim='station', how='all')
                stations = ds.coords["station"].values
                lon_indices = lon_indices.sel(station=stations)
                lat_indices = lat_indices.sel(station=stations)
                stations = metaInfo[metaInfo['id'].isin(stations)]['file_name'].values
                return lat_indices, lon_indices, stations

        def read_netcdfs(soilmositureFiles, lon_indices, lat_indices, step_size=1):
            data = xr.open_mfdataset(soilmositureFiles, decode_cf=True,
                                   concat_dim='time', engine="h5netcdf", 
                                   combine='nested', combine_attrs='drop',
                                   join='override', parallel=True,
                                   mask_and_scale=False, 
                                   chunks={})['upper_soil_saturation_degree']
            # data = data.isel(time=slice(0,365*3))
            
            timeSteps = [i for i in range(0, data.time.size+1, step_size)]
            for timeStep in tqdm(timeSteps, desc="Reading Points", disable=False):
                data_sub = data.isel(time=slice(timeStep,timeStep+step_size)).load().isel(lon=lon_indices, lat=lat_indices)
                if timeStep == 0: dsList = [data_sub]
                else: dsList.append(data_sub)
            selectedData = xr.concat(dsList, 'time')
            return selectedData

        def read_layerDepth(metaInfo, layerDepthFile, lat_indices, lon_indices):
            
            for var in ['firstStorDepth', 'layerDepth']:
                try:
                    data = xr.open_dataset(layerDepthFile, decode_cf=True, engine="h5netcdf", mask_and_scale=False)[var]
                    if var == 'layerDepth':
                        data = data.rename({'latitude': 'lat', 'longitude': 'lon'})
                        data= data.rename('firstStorDepth')
                    break  
                except KeyError:
                    continue  # If the dataset couldn't be opened, continue to the next variable
                
            data = data.isel(lon=lon_indices, lat=lat_indices)
            data = data.to_dataframe().reset_index()
            metaInfo = pd.merge(metaInfo, data, left_on='id', right_on='station')
            metaInfo['depth_end'] = metaInfo['depth_end'].astype(float)
            metaInfo['firstStorDepth'] = metaInfo['firstStorDepth'].astype(float)
            metaInfo = metaInfo[metaInfo['depth_end'] <= metaInfo['firstStorDepth']]
            return metaInfo
        
        search_criteria = create_search_criteria(soilmositureFiles)
        metaInfo = findFluxSites(fluxMetaFile, search_criteria)
        if len(metaInfo) != 0:
            lat_indices, lon_indices, stations = getIndices(soilmositureFiles, metaInfo)
            
            if len(stations) != 0:
                metaInfo = read_layerDepth(metaInfo, layerDepthFile, lat_indices, lon_indices)
                sim_ds = read_netcdfs(soilmositureFiles, lon_indices, lat_indices)
                (saveFolder / "timeseries").mkdir(parents=True, exist_ok=True)	
                (saveFolder / "points").mkdir(parents=True, exist_ok=True)	
                (saveFolder / "scores").mkdir(parents=True, exist_ok=True)	
                
                stationsGrouuped = metaInfo['name'].unique()
                points = pd.DataFrame(columns=['name', 'latitude', 'longitude'])
                for stationGroup in tqdm(stationsGrouuped, desc="Processing stations"):
                    subMeta = metaInfo[metaInfo['name'] == stationGroup]
                    
                    fileNames = subMeta['file_name'].values
                    count = 0
                    for fileName in fileNames:
                        count = count + 1
                        filePath = [s for s in fluxFiles if fileName in str(s)][0]
                        df = pd.read_csv(filePath, delim_whitespace=True, skiprows=1,header=None, names=['date', 'time', 'obs_sm'], usecols=[0,1,2])
                        df['TIMESTAMP'] = pd.to_datetime(df['date'] + ' ' + df['time'])
                        df.set_index('TIMESTAMP', inplace=True)
                        df = df.drop(['date', 'time'], axis=1)
                        frequency = df.index.inferred_freq
                        df = df.resample('D').sum(min_count=1)
                        if count == 1:
                            df_all = df
                        else:
                            df_all = pd.concat([df_all, df], axis=1)
                    df_all = df_all.mean(axis=1)
                    metaInfoSubset = subMeta.iloc[0][['id', 'name', 'latitude', 'longitude']]
                    df_all = df_all.reset_index()
                    df_all['id'] = metaInfoSubset['id']
                    df_all['name'] = metaInfoSubset['name']
                    df_all['latitude'] = metaInfoSubset['latitude']
                    df_all['longitude'] = metaInfoSubset['longitude']
                    df_all = df_all.set_index('TIMESTAMP')
                    df_all = df_all.rename(columns={0: 'obs_sm'})
                    if frequency != 'D' or frequency != 'H':
                        df_all.index = df_all.index.date
                    df_all = pd.DataFrame({'name': df_all['name'].reindex(sim_ds.time.values, fill_value=df_all['name'].values[0]), 
                                        'latitude': df_all['latitude'].reindex(sim_ds.time.values, fill_value=df_all['longitude'].values[0]), 
                                        'longitude': df_all['longitude'].reindex(sim_ds.time.values, fill_value=df_all['latitude'].values[0]),
                                        'obs_sm': df_all['obs_sm'].reindex(sim_ds.time.values),
                                        'sim_sm': sim_ds.sel(station=df_all['id'].values[0]).values})
                    first_non_missing, last_non_missing = df_all['obs_sm'].first_valid_index(), df_all['obs_sm'].last_valid_index()
                    df_all = df_all.loc[first_non_missing:last_non_missing]
                    df_all.to_parquet(saveFolder / f"timeseries/{df_all['name'].values[0]}.parquet")
                    points = pd.concat([points, pd.DataFrame({'name': [metaInfoSubset['name']],
                                                            'latitude': [metaInfoSubset['latitude']],
                                                            'longitude': [metaInfoSubset['longitude']]})], axis=0)

                # Assuming 'points' is your DataFrame
                points['geometry'] = points.apply(lambda row: Point(row.longitude, row.latitude), axis=1)
                geo_points = gpd.GeoDataFrame(points, geometry='geometry')
                geo_points = geo_points.drop(['latitude', 'longitude'], axis=1)
                geo_points.to_file(saveFolder / f"points/points.shp")

    def scoreTimeSeries(saveFolder):
        saveFolder = Path(saveFolder)
        timeSeriesFiles = sorted(saveFolder.glob(f"timeseries/*.parquet"))
        if len( timeSeriesFiles) != 0:
            df = pd.DataFrame(columns=['name', 'KGE', 'KGE_r', 'KGE_alpha', 'KGE_beta'])
            for file in tqdm(timeSeriesFiles, desc="Scoring Stations Daily", disable=False):
                file = pd.read_parquet(file)
                station=file['name'].values[0]
                kge, r, alpha, beta = he.evaluator(he.kge, file["sim_sm"], file["obs_sm"])
                if df.empty:
                    df = pd.DataFrame({'name': [station], 
                                       'KGE': kge,
                                       'KGE_r': r,
                                       'KGE_alpha': alpha,
                                       'KGE_beta': beta})
                else:
                    df = pd.concat([df, pd.DataFrame({'name': [station], 
                                                      'KGE': kge,
                                                      'KGE_r': r,
                                                      'KGE_alpha': alpha,
                                                      'KGE_beta': beta})], ignore_index=True)
            scores = gpd.read_file(saveFolder / "points/points.shp")
            scores = scores.merge(df, on="name")
            scores.to_file(saveFolder / f"scores/day_scores.shp")
            
            df = pd.DataFrame(columns=['name', 'KGE', 'KGE_r', 'KGE_alpha', 'KGE_beta'])
            for file in tqdm(timeSeriesFiles, desc="Scoring Stations Monthly", disable=False):
                file = pd.read_parquet(file)
                station=file['name'].values[0]
                file = file.resample('M').agg({'name': 'first',
                                                'latitude': 'first',
                                                'longitude': 'first',
                                                'obs_sm': 'mean',
                                                'sim_sm': 'mean'})
                kge, r, alpha, beta = he.evaluator(he.kge, file["sim_sm"], file["obs_sm"])
                if df.empty:
                    df = pd.DataFrame({'name': [station], 
                                       'KGE': kge,
                                       'KGE_r': r,
                                       'KGE_alpha': alpha,
                                       'KGE_beta': beta})
                else:
                    df = pd.concat([df, pd.DataFrame({'name': [station], 
                                                      'KGE': kge,
                                                      'KGE_r': r,
                                                      'KGE_alpha': alpha,
                                                      'KGE_beta': beta})], ignore_index=True)
            scores = gpd.read_file(saveFolder / "points/points.shp")
            scores = scores.merge(df, on="name")
            scores.to_file(saveFolder / "scores/month_scores.shp")
            print('done')