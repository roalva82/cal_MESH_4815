#!/usr/bin/env python
# coding: utf-8

import subprocess
import os
import pandas as pd
import numpy as np
import csv
import xarray as xr

os.chdir('./model')

def nse(simulated, observed):
	denominator = np.sum((observed - np.mean(observed)) ** 2)
	numerator = np.sum((simulated - observed) ** 2)
	nse = 1 - numerator / denominator
	if np.isnan(nse):
		#pass a nse of 0 to avoid adding it to the rest and keeps the counter at 0
		return 0,0
	else:
		return nse,1

#drainage database file
drainage = xr.open_dataset('./MESH_drainage_database.nc')

#read list of stations with data
stations = pd.read_csv('./obsFiles/stationsMapping.csv')

#path to observed files
path = './obsFiles/hydrometric-daily-mean/'

#range of dates to compute metric
start_date = '2015-03-01'
end_date = '2015-03-31'

# run mesh
mesh_command  = './sa_mesh'
subprocess.run([mesh_command])

#add rank as a column to the stations mappings dataframe
stations['rank'] = ""

#map the rank to the specific station
for el in range(len(stations)):
	sb = stations['COMID'][el]
	stations.loc[el,'rank'] = drainage.sel(subbasin=sb)['Rank'].values

#read simulated data
sim = pd.read_csv('./results/ROF_D.csv', header=None)
sim = sim.rename(columns={0: 'DATE'})
sim['DATE'] = pd.to_datetime(sim['DATE'])

metric_total = 0
valid_obs = 0

#loop over gauging stations and obtain averaged metric
for el in range(len(stations)):
	#extract the simulated time series
	mask_sim = (sim['DATE'] >= start_date) & (sim['DATE'] <= end_date)
	rank = stations.loc[el,'rank'].item()
	sim_masked = sim.loc[mask_sim][['DATE',rank]]

	#read the observed time series
	filename = path + stations.loc[el,'ID'] + '_DISCHARGE.csv'
	obs = pd.read_csv(filename)

	#extract data for the range of dates
	obs['DATE'] = pd.to_datetime(obs['DATE'])
	mask_obs = (obs['DATE'] >= start_date) & (obs['DATE'] <= end_date)
	obs_masked = obs.loc[mask_obs][['DATE','DISCHARGE']]

	#merge for values that exist on both time series
	merged_df = pd.merge(sim_masked, obs_masked, on='DATE', how='inner')

	#extract simulated and observed timeseries
	sim_trim = merged_df.iloc[:,1].to_numpy()
	obs_trim = merged_df.iloc[:,2].to_numpy()

	if sim_trim.size != 0 or obs_trim.size != 0:
		metric, counter = nse(sim_trim,obs_trim)
		metric_total += metric
		valid_obs += counter

#compute average of metric (multiply by -1 for nse)
final_metric = -1.0*metric_total/valid_obs

#write results to metric file
with open('./results/Metric.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow([final_metric])
    

