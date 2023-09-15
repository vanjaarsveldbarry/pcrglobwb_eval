Inspired by and copied from https://github.com/JannisHoch/pcrglobwb_utils

## Validate Daily Discharge
1.) Download Daily Discharge Data from the [GRDC Data Portal](https://portal.grdc.bafg.de/applications/public.html?publicuser=PublicUser#dataDownload/Stations)

2.) Run grdcPre.py once to process grdc format into usable format

2.) Run grdcEval.py to get:

*Validation Dataset*
```python
<xarray.Dataset>
Dimensions:              (station: 296, time: 5)
Coordinates:
    cat_area             (station) float64 1.561e+03 349.8 ... 66.5 73.3
    lat                  (station) float64 49.76 50.04 48.8 ... 46.49 46.57
    lon                  (station) float64 16.98 16.91 16.85 ... 9.898 9.936
    name                 (station) <U28 'MORAVICANY (No. 6142100)' ... 'LA PU...
  * station              (station) int64 6142100 6142101 ... 6943160 6943170
  * time                 (time) datetime64[ns] 1979-01-01 ... 1979-01-05
    ts_end               (station) datetime64[ns] 2022-12-01 ... 2019-01-01
    ts_start             (station) datetime64[ns] 1911-11-01 ... 1954-05-01
Data variables:
    grdc_discharge       (station, time) float64 76.6 55.7 50.3 ... 0.64 0.639
    mean_grdc_discharge  (station) float64 17.2 6.093 33.39 ... 2.821 2.265
    sim_discharge        (station, time) float32 0.01895 2.133 ... 0.6857 0.67
```

*Scores Dataset*
```python
<xarray.Dataset>
Dimensions:           (station: 296)
Coordinates:
    cat_area          (station) float64 1.561e+03 349.8 1.228e+04 ... 66.5 73.3
    lat               (station) float64 49.76 50.04 48.8 ... 46.49 46.49 46.57
    lon               (station) float64 16.98 16.91 16.85 ... 9.906 9.898 9.936
    name              (station) <U28 'MORAVICANY (No. 6142100)' ... 'LA PUNT ...
  * station           (station) int64 6142100 6142101 ... 6943160 6943170
    ts_end            (station) datetime64[ns] 2022-12-01 ... 2019-01-01
    ts_start          (station) datetime64[ns] 1911-11-01 ... 1954-05-01
Data variables:
    kge               (station) float64 -0.6826 -0.7006 nan ... -312.9 -0.831
    kge_pearson_coef  (station) float64 -0.5011 -0.5273 nan ... -0.1695 0.9695
    kge_rel_var       (station) float64 1.359 0.8259 nan ... 629.1 314.5 2.829
    kge_bias          (station) float64 0.3299 0.2726 inf ... 8.224 15.88 1.085
    rsquared          (station) float64 0.2511 0.2781 nan ... 0.02872 0.9399
    mse               (station) float64 9.842e+03 1.207e+03 ... 31.95 0.01654
    rmse              (station) float64 99.21 34.74 0.0 ... 5.343 5.652 0.1286
```
