# Validation Metrics
# JCS
# October 3, 2023

from dataclasses import dataclass

from calendar import monthrange 
import numpy as np
import warnings
from scipy.stats.stats import spearmanr
from scipy.stats import cumfreq
from mpl_toolkits.axes_grid1 import make_axes_locatable
import xarray as xr

# Metrics we are using:

def bias(obs, mod):
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  out = np.mean(obs[sel]) -np.mean(mod[sel])
  return out, len(sel)

def kge(obs, mod):
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  cc = spearmanr(obs[sel],mod[sel])[0,1]
  alpha = np.std(obs[sel])/np.std(mod[sel])
  beta  = np.sum(obs[sel])/np.sum(mod[sel])
  kge   = 1- np.sqrt( (cc-1)**2 + (alpha-1)**2 + (beta-1)**2 )
  return kge, cc, alpha, beta

def normalizeMonth(data):
  seasonCycle = np.zeros(len(data))
  for m in range(12):
    monthSelection = np.arange(m,len(data), 12)
    seasonCycle[monthSelection] = np.mean(data[monthSelection])
  return data - seasonCycle

def anomalyCorrelation(obs, mod, timeScale = "month"):
  normObs = normalizeMonth(obs)
  normMod = normalizeMonth(mod)
  return spearmanr(normObs, normMod)[0]







### Other metrics!
# def nashSutcliffe(obs, mod):
#   obsSel = np.isnan(obs) == False
#   modSel = np.isnan(mod) == False
#   sel = obsSel & modSel
#   MSE = np.sum((obs[sel]-mod[sel])**2)
#   MSEclim = np.sum((obs[sel] - np.mean(obs[sel]))**2)
#   NSE = 1-MSE/MSEclim
#   return np.maximum(NSE,-100.)

# def rmse(obs, mod):
#   obsSel = np.isnan(obs) == False
#   modSel = np.isnan(mod) == False
#   sel = obsSel & modSel
#   out = np.mean((obs[sel]-mod[sel])**2)**0.5
#   return out

# def calculateMetrics(obs, mod):
#   timeSize = 12
#   obsSel = np.isnan(obs) == False
#   modSel = np.isnan(mod) == False
#   sel = obsSel & modSel
#   #print(obs[sel])
#   #print(mod[sel])
#   if len(obs[sel]) > timeSize and len(mod[sel]) > timeSize:
#     R = spearmanr(obs[sel], mod[sel])[0]
#     NS = nashSutcliffe(obs[sel], mod[sel])
#     RMSE = rmse(obs[sel], mod[sel])
#     Bias, numPoints = bias(obs[sel], mod[sel])
#     KGE, CC, Alpha, Beta = kge(obs[sel], mod[sel])
#     AC = anomalyCorrelation(obs[sel], mod[sel])
#     #AC =0
#     #Bias, numPoints = 0
#     #RMSE = 0
#     #NS =0
#     return R, AC, KGE, CC, Alpha, Beta, NS, RMSE, Bias, numPoints
#   else:
#     return np.zeros((10))
