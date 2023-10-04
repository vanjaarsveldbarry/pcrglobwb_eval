from dataclasses import dataclass
from pathlib import Path
import xarray as xr
# import geopandas as gpd
# from rasterio import features
import numpy as np
# import hydromt as hmt
import xskillscore as xs


@dataclass        
class snow:
    """Working with grdc data.
    """
   
    snowDataDirectory: str
    simDirectory: str
    validationData: xr.Dataset = None
    scores: xr.Dataset = None
    
    
    def __post_init__(self):
        self.snowDataDirectory = Path(self.snowDataDirectory)
        self.simDirectory = Path(self.simDirectory)
        self.validationData = dailyValidate(self.snowDataDirectory,
                                            self.simDirectory)
        self.scores = scores(self.validationData.sim_snow_cover, 
                                   self.validationData.obs_snow_cover)

def scores(df_sim, df_obs):
    dichotomous_category_edges = np.array([0, 0.5, 1])
    dichotomous_contingency = xs.Contingency(
    df_obs, df_sim, dichotomous_category_edges, dichotomous_category_edges, dim=["lat", "lon", "time"])
    print(dichotomous_contingency.bias_score())
    
    
    return scores_ds
    
def dailyValidate(snowObservedDataDirectory, simDirectory):
    """Array-wise validation.
    """
    
    def readSim(simDirectory):
        sweFiles = sorted(simDirectory.glob(f"**/*netcdf/snowCoverSWE_dailyTot_output.nc"))
        
        pcr_ds = xr.open_mfdataset(sweFiles, engine='netcdf4')
        #HACK - must remove
        pcr_ds = pcr_ds.sel(time=slice("2000-02-24", "2000-02-29")).compute()
        #######
        pcr_ds = xr.where(pcr_ds > 0, 1., pcr_ds.snow_water_equivalent)   
        return pcr_ds.rename({'snow_water_equivalent': 'sim_snow_cover'})
    
    def readObs(obsFile, sim_ds):
        with xr.open_zarr(obsFile / 'snowCover.zarr', chunks='auto') as snowData:
            lonMin, lonMax = sim_ds.lon.min().values, sim_ds.lon.max().values
            latMin, latMax = sim_ds.lat.min().values, sim_ds.lat.max().values
            snowData = snowData.reindex_like(sim_ds)
            # snowData = snowData.sel(lon=slice(lonMin.item(), lonMax.item()), 
            #                         lat=slice(latMax.item(), latMin.item()),
            #                         time=slice(sim_ds.time.values[0], sim_ds.time.values[-1]))
        return snowData.rename({'scfg': 'obs_snow_cover'}).compute()
    
    sim_ds = readSim(simDirectory)
    obs_ds = readObs(snowObservedDataDirectory, sim_ds)    
    validation_ds = xr.merge([sim_ds, obs_ds])
    return validation_ds