
class preprocess:
    """Preprocessig input data for further validation.
    """
    def grdcData(self, grdcRawDataDirectory: str, saveFolder: str):
        """Convert grdc data to something useful.
        """
        saveFolder.mkdir(exist_ok=True)
        (saveFolder / "temp").mkdir(exist_ok=True)

        dayFiles = sorted(grdcRawDataDirectory.glob("**/*.txt"))
        totalLength = len(dayFiles)

        count = 0
        for file in tqdm(dayFiles):
            count = count + 1
            grdc_data = pcrglobwb_utils.obs_data.grdc_data(file)
            stationInfo = grdc_data.get_grdc_station_properties()
            data = grdc_data.get_grdc_station_values(var_name="OBS", col_name=" Value")
            data = data["OBS"]

            timeVals = data.index.values
            data = data.values
            data = data.reshape(data.shape + (1,))

            latVals = np.array([stationInfo["latitude"]])
            lonVals = np.array([stationInfo["longitude"]])
            stationVals = np.array([stationInfo["grdc_no"]])
            stationgGrid = np.array([stationInfo["grdc_no"]])
            stationgGrid = stationgGrid.reshape((1,) + (1,) + stationgGrid.shape)

            dataSet = xr.Dataset(
                data_vars={
                    "grdc_discharge": ((["time", "station"]), data),
                },
                coords={
                    "time": (["time"], timeVals),
                    "station": (["station"], stationVals),
                    "lat": (["station"], np.array([stationInfo["latitude"]])),
                    "lon": (["station"], np.array([stationInfo["longitude"]])),
                    "ts_start": (["station"], np.array([stationInfo["ts_start"]])),
                    "ts_end": (["station"], np.array([stationInfo["ts_end"]])),
                    "cat_area": (["station"], np.array([stationInfo["cat_area"]])),
                    "name": (
                        ["station"],
                        np.array([f"{stationInfo['station']} (No. {stationInfo['grdc_no']})"]),
                    ),
                },
            )
            if count == 1:
                full_ds = dataSet

            else:
                dataSet, full_ds = xr.align(
                    dataSet, full_ds, exclude=["lat", "lon"], join="outer"
                )

                full_ds = xr.merge([dataSet, full_ds], join="outer", compat="no_conflicts")

                if not count % 10:
                    full_ds = full_ds.isel(station=slice(None, -1))
                    full_ds.to_zarr(saveFolder / f"temp/test{count}.zarr", mode="w")
                    full_ds = full_ds.isel(station=slice(-2, -1))

        fullFiles = sorted((saveFolder / f"temp").glob("*/"))
        full_ds = xr.open_mfdataset(
            fullFiles, engine="zarr", combine="nested", concat_dim="station"
        )
        full_ds["name"] = full_ds["name"].astype("str")
        full_ds = full_ds.compute()
        full_ds["mean_grdc_discharge"] = full_ds["grdc_discharge"].mean("time")
        full_ds.to_zarr(saveFolder / "GRDC_array.zarr", mode="w", consolidated=True)

        count = 0
        for file in tqdm(dayFiles):
            count = count + 1
            grdc_data = pcrglobwb_utils.obs_data.grdc_data(file)
            stationInfo = grdc_data.get_grdc_station_properties()
            stationgGrid = np.array([stationInfo["grdc_no"]])
            stationgGrid = stationgGrid.reshape((1,) + stationgGrid.shape)

            dataSet = pd.DataFrame(
                {
                    "stationID": np.array([stationInfo["grdc_no"]]),
                    "lat": np.array([stationInfo["latitude"]]),
                    "lon": np.array([stationInfo["longitude"]]),
                }
            )

            if count == 1:
                full_ds = dataSet

            else:
                full_ds = pd.concat([dataSet, full_ds], axis=0, ignore_index=False)

        (saveFolder / "GRDC_points").mkdir(exist_ok=True)
        gdf = gpd.GeoDataFrame(
            full_ds, geometry=gpd.points_from_xy(full_ds.lon, full_ds.lat), crs="EPSG:4326"
        )
        gdf.to_file(saveFolder / "GRDC_points/GRDC_points.shp")
        shutil.rmtree(saveFolder / "temp")
    