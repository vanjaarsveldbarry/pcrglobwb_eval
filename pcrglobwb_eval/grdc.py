from dataclasses import dataclass
from pathlib import Path
import xarray as xr
import geopandas as gpd
from rasterio import features
import numpy as np
import hydromt as hmt
import Other_Metrics

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
    df_sim = df_sim.compute()
    kge_ds = hmt.stats.skills.kge(df_sim, df_obs)
    r2_ds = hmt.stats.skills.rsquared(df_sim, df_obs)
    mse_ds = hmt.stats.skills.mse(df_sim, df_obs)
    rmse_ds = hmt.stats.skills.rmse(df_sim, df_obs)
    non_parametric_kge = Other_Metrics.kge(df_obs, df_sim)
    acc = Other_Metrics.anomalyCorrelation(df_obs, df_sim)
    bias = Other_Metrics.bias(df_obs, df_sim)
    
    scores_ds = xr.merge([kge_ds, r2_ds, mse_ds, rmse_ds, acc, non_parametric_kge, bias])
    return scores_ds
    
def dailyValidate(grdcObservedDataDirectory, simDirectory):
    """Point-wise validation.
    """
    def readSim(simDirectory):    
        dischargeFiles = sorted(simDirectory.glob(f"**/*discharge_dailyTot_output.nc"))
        pcr_ds = xr.open_mfdataset(dischargeFiles, engine='netcdf4')
        #HACK - must remove
        pcr_ds = pcr_ds.sel(time=slice("1995-01-01", "1995-01-05")).compute()
        pcr_ds["mean_discharge"] = pcr_ds.discharge.mean("time")
        return pcr_ds
    
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
            return gdf
        
        # Create mask and select GRDC points
        mask = sim_ds.isel(time=0).discharge
        mask = mask.rio.set_crs("epsg:4326")
        transform = mask.rio.transform()
        mask = xr.where(mask >= 0., 1, 0)
        mask = mask.astype('uint8').values
        mask = vectorize(mask, 0, transform, crs="epsg:4326", name="value")
        mask = mask.dissolve(by='value')
            
        grdcPoints = gpd.read_file(obsFile / 'GRDC_points/GRDC_points.shp', engine='pyogrio')
        grdcPoints = grdcPoints.overlay(mask, how='intersection')
        grdcPoints = grdcPoints.stationID.unique()

        #clip GRDC data to selected GRDC points
        with xr.open_zarr(obsFile / 'GRDC_array.zarr', chunks='auto') as grdcData:
            grdcData = grdcData.drop_duplicates(dim="station")
            grdcData = grdcData.sel(station=grdcData.station.isin(grdcPoints), 
                                    time=slice(sim_ds.time.values[0], sim_ds.time.values[-1]))
        return grdcData.compute()
    
    def add_sim_discharge(ds):
        window = 5
        tolerance = window * 0.008333333
        lonMin, lonMax = ds.lon.values - tolerance, ds.lon.values + tolerance
        latMin, latMax = ds.lat.values - tolerance, ds.lat.values + tolerance
        sub_pcr = sim_ds.sel(lon=slice(lonMin.item(), lonMax.item()), lat=slice(latMax.item(), latMin.item()))
        sub_ID = sub_pcr.mean_discharge - ds.mean_grdc_discharge.values
        sub_ID = np.absolute(sub_ID)
        xyID = sub_ID.argmin(dim=["lat", "lon"])
        y, x = np.array([xyID["lat"].values]), np.array([xyID["lon"].values])
        sub_pcr = sub_pcr.discharge.isel(lat=y, lon=x).squeeze(["lat", "lon"])
        ds["sim_discharge"] = sub_pcr
        return ds
    
    sim_ds = readSim(simDirectory)
    obs_ds = readObs(grdcObservedDataDirectory, sim_ds)
    validation_ds = obs_ds.groupby("station").map(add_sim_discharge)
    return validation_ds