from dataclasses import dataclass
from pathlib import Path
import xarray as xr
import geopandas as gpd
from rasterio import features
import numpy as np
import hydromt as hmt
import warnings

@dataclass
class grdc:
    """Working with gwdepth data.
    """
    gwDepthDataDirectory: str
    simDirectory: str
    validationData: xr.Dataset = None
    scores: xr.Dataset = None

    def __post_init__(self):
        self.gwDepthDataDirectory = Path(self.gwDepthDataDirectory)
        self.simDirectory = Path(self.simDirectory)
        self.validationData = self.monthly_validate()
        self.scores = self.scores(self.validationData.gwDepth, 
                                  self.validationData.sim_ds)

    def scores(df_sim, df_obs):
        nse_ds = hmt.stats.skills.nashsutcliffe(df_sim, df_obs)
        kge_ds = hmt.stats.skills.kge(df_sim, df_obs)
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


    def monthly_validate(self):
        def extract_sim_data():
            gwDepthFiles = sorted(self.simDirectory.glob("**/*groundwaterDepthLayer1_monthEnd_output.nc"))
            pcr_gw = xr.open_mfdataset(gwDepthFiles, combine ="by_coords",engine='netcdf4', chunks=None)     
            return pcr_gw.compute()

        def read_obs(obs_file, sim_gw):
            def vectorize(data, nodata, transform, crs="epsg:4326", name="value"):
                feats_gen = features.shapes(data, mask=data != nodata, transform=transform, connectivity=8)
                feats = [{"geometry": geom, "properties": {name: val}} for geom, val in list(feats_gen)]

                # # parse to geopandas for plotting / writing to file
                gdf = gpd.GeoDataFrame.from_features(feats, crs=crs)
                gdf[name] = gdf[name].astype(data.dtype)
                return gdf

            mask = sim_gw.isel(time=0).groundwater_depth_for_layer_1
            mask = mask.rio.set_crs("epsg:4326")
            transform = mask.rio.transform()
            mask = xr.where(mask >= 0., 1, 0)
            mask = mask.astype('uint8').values
            mask = vectorize(mask, 0, transform, crs="epsg:4326", name="value")
            mask = mask.dissolve(by='value')

            gw_depth_points = gpd.read_file(obs_file / 'well_points/obswell_points.shp', engine='pyogrio')
            gw_depth_points = gw_depth_points.overlay(mask, how='intersection')
            gw_depth_points = gw_depth_points.well.unique()
            # with xr.open_zarr(obs_file / 'gwDepth_array.zarr', chunks='auto') as obs_well_data:
            #         obs_well_data = obs_well_data.sel(
            #             well=obs_well_data.well.isin(gw_depth_points),
            #             time=slice(sim_gw.time.values[0], sim_gw.time.values[-1])
            #     )
            #         print(obs_well_data)
            # return obs_well_data.compute()


            # Open your Zarr dataset
            with xr.open_zarr(obs_file / 'gwDepth_array.zarr', chunks='auto') as obs_well_data:
                # Define the time range you want to select
                start_time = '2000-01-01'  # Replace with the start date in ISO format
                end_time = '2016-12-31'    # Replace with the end date in ISO format

                # Use .sel() to select the time frame
                obs_well_data = obs_well_data.sel(
                    well=obs_well_data.well.isin(gw_depth_points),
                    time=slice(start_time, end_time)
                )

            return obs_well_data
            
        def add_sim_gw_depth(gw):
            # Extract the latitude and longitude coordinates from sub_pcr
            lat_coords = gw.lat.values
            lon_coords = gw.lon.values

            # Use these coordinates to select values from sim_gw
            sub_pcr = sim_gw.sel(lat=lat_coords, lon=lon_coords, method="nearest")
            
            # Assuming you want to extract "groundwater_depth_for_layer_1" from sub_pcr
            sub_pcr = sub_pcr.groundwater_depth_for_layer_1.squeeze()

            # Assign the selected data to "sim_groundwater_depth1" in gw
            gw["sim_ds"] = sub_pcr
            

            return gw.compute()


        sim_gw = extract_sim_data()
        obs_gw = read_obs(self.gwDepthDataDirectory, sim_gw)
        validation_gw = obs_gw.groupby("well").map(add_sim_gw_depth)
        print(validation_gw)   
        return validation_gw
 

    
    
    
    
    
    
    
    
    
    
    
