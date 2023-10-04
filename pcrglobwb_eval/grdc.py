from dataclasses import dataclass
from pathlib import Path
import xarray as xr
import geopandas as gpd
from rasterio import features
import numpy as np
import hydromt as hmt

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

@dataclass        
class grdc:
    """Working with grdc data.
    """
   
    grdcDataDirectory: str
    simDirectory: str
    validationData: xr.Dataset = None
    scores: xr.Dataset = None
    
    
    def __post_init__(self):
        self.grdcDataDirectory = Path(self.grdcDataDirectory)
        self.simDirectory = Path(self.simDirectory)
        
        self.validationData = dailyValidate(self.grdcDataDirectory,
                                            self.simDirectory)
        
        self.scores = scores(self.validationData.sim_discharge, 
                             self.validationData.grdc_discharge)

def scores(df_sim, df_obs):
    nse_ds = hmt.stats.skills.nashsutcliffe(df_sim, df_obs)
    kge_ds = hmt.stats.skills.kge(df_sim, df_obs)
    #TODO non parametric kge not working
    # kge_np_ds = hmt.stats.skills.kge_non_parametric(df_sim, df_obs)
    bias_ds = hmt.stats.skills.bias(df_sim, df_obs)
    bias_percent_ds = hmt.stats.skills.percentual_bias(df_sim, df_obs)
    spearman_ds = hmt.stats.skills.spearman_rank_correlation(df_sim, df_obs)
    r2_ds = hmt.stats.skills.rsquared(df_sim, df_obs)
    mse_ds = hmt.stats.skills.mse(df_sim, df_obs)
    rmse_ds = hmt.stats.skills.rmse(df_sim, df_obs)
    def acc_calc(df_sim, df_obs):
        sim_mean = df_sim.mean(dim='time')
        obs_mean = df_obs.mean(dim='time')
        simAnom = df_sim - sim_mean
        obsAnom = df_obs - obs_mean
        acc_ds = hmt.stats.skills.spearman_rank_correlation(df_sim, df_obs) 
        acc_ds = acc_ds.rename('acc_coeff')
        return acc_ds
    acc_ds = acc_calc(df_sim, df_obs)
    scores_ds = xr.merge([nse_ds,kge_ds,bias_ds,bias_percent_ds,
                          spearman_ds,r2_ds, mse_ds, rmse_ds, acc_ds])
    return scores_ds
    
def dailyValidate(grdcObservedDataDirectory, simDirectory):
    """Point-wise validation.
    """
    def readSim(simDirectory):    
        dischargeFiles = sorted(simDirectory.glob(f"**/*netcdf/discharge_dailyTot_output.nc"))
        pcr_ds = xr.open_mfdataset(dischargeFiles, engine='h5netcdf', chunks=None)
        #HACK - must remove
        pcr_ds = pcr_ds.sel(time=slice("1979-01-01", "1979-01-20"))
        ######
        # pcr_ds["mean_discharge"] = pcr_ds.discharge.mean("time")
        return pcr_ds.compute()
    
    def readUpstreamArea(obsFile, sim_ds):
        with xr.open_zarr(obsFile / 'GRDC_upstreamArea.zarr', chunks='auto') as up_ds:
            up_ds = up_ds.sel(longitude=slice(sim_ds.lon.min(), sim_ds.lon.max()),
                              latitude=slice(sim_ds.lat.max(), sim_ds.lat.min()))
        return up_ds.compute()
    
    def readObs(obsFile, sim_ds):
        
        def vectorize(data, nodata, transform, crs="epsg:4326", name="value"):
            feats_gen = features.shapes(data,
                                        mask=data != nodata,
                                        transform=transform,
                                        connectivity=8)
            feats = [{"geometry": geom, "properties": {name: val}} for geom, val in list(feats_gen)]

            # # parse to geopandas for plotting / writing to file
            gdf = gpd.GeoDataFrame.from_features(feats, crs=crs)
            gdf[name] = gdf[name].astype(data.dtype)
            return gdf.dissolve(by='value')
        
        # Create mask and select GRDC points
        mask = sim_ds.isel(time=0).discharge
        mask = xr.where(mask >= 0., 1, 0)
        mask = mask.rio.set_crs("epsg:4326")
        transform = mask.rio.transform()
        mask = mask.astype('uint8').values
        mask = vectorize(mask, 0, transform, crs="epsg:4326", name="value")
        grdcPoints = gpd.read_file(obsFile / 'GRDC_points/GRDC_points.shp', engine='pyogrio')
        grdcPoints = grdcPoints.overlay(mask, how='intersection')
        grdcPoints = grdcPoints.stationID.unique()

        #clip GRDC data to selected GRDC points
        with xr.open_zarr(obsFile / 'GRDC_array.zarr', chunks='auto') as grdcData:
            #HACK checck that the processed data npo longer contains duplicate station ID and \
            #     figure out how to handle them 
            grdcData = grdcData.drop_duplicates(dim="station")
            grdcData = grdcData.sel(station=grdcData.station.isin(grdcPoints), 
                                    time=slice(sim_ds.time.values[0], sim_ds.time.values[-1]))
            
        return grdcData.compute()
    
    def add_sim_discharge(ds):
        #5 km window
        window = 5 * 0.008333333
        lonMin, lonMax = (ds.lon.values - window).item(), (ds.lon.values + window).item()
        latMin, latMax = (ds.lat.values - window).item(), (ds.lat.values + window).item()
        sub_upstream_area = upstreamArea.sel(longitude=slice(lonMin, lonMax), 
                                             latitude=slice(latMax, latMin))
        diff = np.absolute(sub_upstream_area.upstream_area - ds.cat_area.values.item())
        xyID = diff.argmin(dim=["latitude", "longitude"])
        sub_upstream_area = sub_upstream_area.upstream_area.isel(latitude=np.array([xyID["latitude"].values]), 
                                                                 longitude=np.array([xyID["longitude"].values])).squeeze(["latitude", "longitude"])
        sub_pcr = sim_ds.discharge.sel(lat=sub_upstream_area.latitude.values.item(), 
                                        lon=sub_upstream_area.longitude.values.item(), method = 'nearest')#.squeeze(["latitude", "longitude"])
        ds["sim_discharge"] = sub_pcr
        # ds["mean_sim_discharge"] = sub_pcr.mean(dim='time')
        return ds

    sim_ds = readSim(simDirectory)
    upstreamArea = readUpstreamArea(grdcObservedDataDirectory, sim_ds)
    obs_ds = readObs(grdcObservedDataDirectory, sim_ds)
    obs_ds = obs_ds.isel(station=slice(0,1))
    start_timeF = time.time()
    validation_ds = obs_ds.groupby("station").map(add_sim_discharge, shortcut =True)
    end_timeF = time.time()
    execution_timeF = end_timeF - start_timeF
    print(f"WHOLE took {execution_timeF:.6f} seconds to run")
    return validation_ds.compute()