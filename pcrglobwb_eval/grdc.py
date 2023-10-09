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

    def scores(self, df_sim, df_obs):
        df_sim = df_sim.compute()
        kge_gw = hmt.stats.skills.kge(df_sim, df_obs)
        r2_gw = hmt.stats.skills.rsquared(df_sim, df_obs)
        mse_gw = hmt.stats.skills.mse(df_sim, df_obs)
        rmse_gw = hmt.stats.skills.rmse(df_sim, df_obs)
        scores_gw = xr.merge([kge_gw, r2_gw, mse_gw, rmse_gw])
        return scores_gw

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

            try:
                with xr.open_zarr(obs_file / 'gwDepth_array.zarr', chunks='auto') as obs_well_data:
                    obs_well_data = obs_well_data.sel(
                        well=obs_well_data.well.isin(gw_depth_points),
                        time=slice(sim_gw.time.values[0], sim_gw.time.values[-1])
                )
                print(obs_well_data)
                return obs_well_data.compute()
            except Exception as e:
            # If an exception is raised (e.g., if the dataset isn't successfully created), show a warning
                warning_message = f"An error occurred when reading observed data: {str(e)}"
                warnings.warn(warning_message)
                return None 

        def add_sim_gw_depth(gw):
            window = 1
            tolerance = window * 0.008333333
            lon_min, lon_max = gw.lon.values - tolerance, gw.lon.values + tolerance
            lat_min, lat_max = gw.lat.values - tolerance, gw.lat.values + tolerance
            sub_pcr = sim_gw.sel(lon=slice(lon_min.item(), lon_max.item()), lat=slice(lat_max.item(), lat_min.item()))
            a = sub_pcr.groundwater_depth_for_layer_1
            a = a.values.reshape(-1, 2)
            b = gw.gwDepth.values
            print(b)
            sub_ID = a - gw.gwDepth.values
            sub_ID = np.absolute(sub_ID)
            xy_ID = sub_ID.argmin(dim=["lat", "lon"])
            y, x = np.array([xy_ID["lat"].values]), np.array([xy_ID["lon"].values])
            sub_pcr = sub_pcr.groundwater_depth1.isel(lat=y, lon=x).squeeze(["lat", "lon"])
            gw["sim_groundwater_depth1"] = sub_pcr
            return gw

        sim_gw = extract_sim_data()
        obs_gw = read_obs(self.gwDepthDataDirectory, sim_gw)
        validation_gw = obs_gw.groupby("well").map(add_sim_gw_depth)
        return validation_gw
