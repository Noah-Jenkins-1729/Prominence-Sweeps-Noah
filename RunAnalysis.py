# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:41:55 2024

@author: Hannah
"""
import sys
import numpy as np
# append necessary file paths, and change E -> D or vice versa
from MeasurementInfo import MeasurementInfo
from RunInfo import RunInfo
import heapq
from scipy import signal
from scipy.optimize import curve_fit
import AnalyzePDE
from AnalyzePDE import SPE_data
from AnalyzePDE import Alpha_data
import matplotlib.pyplot as plt
import matplotlib as mpl
import ProcessWaveforms_MultiGaussian
from ProcessWaveforms_MultiGaussian import WaveformProcessor as WaveformProcessor
import csv
import os
import pandas as pd
#1us, 10x gain, filter on - see Noah's report on OneDrive
invC_filter = 0.0115
invC_err_filter = 0.000098



df_redchi = pd.read_csv('20250715 Prominence Sweep; Data = October 17, 2024 LED SPE; redchi.txt',index_col=0)
PINTERCEPT = np.array(df_redchi.index.values)
PSLOPE = np.array(df_redchi.columns.astype('float').values)

n=30
m=0
bprom = PINTERCEPT[n] #intercept at bias = 32
mprom = PSLOPE[m]  #slope of prominence line


# loop quickly gets the file names
# separate files for each bias voltage -> specifyAcquisition = False
path ='20241017_SPE_LED_405_LXe/'
files = os.listdir(path)
# Filtering only the files
files = [f for f in files if os.path.isfile(path+'/'+f)]
print(files)
#%% read minimal waveforms as a shortcut to print all the metadata
runs = []
for file in range(len(files)):
    run_spe = RunInfo([path+files[file]], specifyAcquisition = False, do_filter = True, upper_limit = 5, poly_correct=True, baseline_correct = True, prominence = 0.008, is_led = True, num_waveforms = 1)
    runs.append(run_spe)
biases = [run.bias for run in runs]
print(biases)
scoperange = [run.yrange for run in runs]
print(scoperange)
offset = [run.offset for run in runs]
print(offset)
#%%
# in general upper_limit is best set at scoperange - offset - 0.001; not always
upperlim = np.array(scoperange)-np.array(offset)-0.002
#%%  
upperlim[0]=0.119 #wasn't getting rid of the scintillation so had to set upper_limit lower
#upperlim = upperlim[:-1]
runs = []
#print(files)
#print(biases)
for file in range(len(files)): #skipping the last run because it happens to be the noise pedestal in this case
    prom = mprom * ( biases[file] - 32 ) + bprom
    run_spe = RunInfo([path+files[file]], specifyAcquisition = False, do_filter = True, upper_limit = upperlim[file], poly_correct=True, baseline_correct = True, prominence = prom, is_led = True)
    if run_spe.bias >30:
        runs.append(run_spe)
#%% baseline noise - 20V
f = '20241017_SPE_LED_405_LXe/Run_1729184064.hdf5'
pedestal = RunInfo([f],  specifyAcquisition = False, do_filter = True, is_solicit = True, upper_limit = 0.008, baseline_correct = True,  poly_correct=False)
# this cell defines the function to get the subtracted histogram
def get_subtract_hist_mean(bias, data1, data2, numbins = 400, plot = False):
    if plot:
        plt.figure()
        (n, b, p) = plt.hist(data1, bins = numbins, density = False, label = 'LED-On', histtype='step')
        plt.axvline(x = np.mean(data1), color = 'blue')
        print('LED on hist: ' + str(np.mean(data1)))
        print('LED off hist: ' + str(np.mean(data2)))
        plt.axvline(x = np.mean(data2), color = 'orange')
        plt.hist(data2, bins = b, density = False, label = 'LED-Off', histtype='step')
    counts1, bins1 = np.histogram(data1, bins = numbins, density = False)
    counts2, bins2 = np.histogram(data2, bins = bins1, density = False)
    centers = (bins1[1:] + bins1[:-1])/2
    subtracted_counts = counts1 - counts2
    # subtracted_counts[subtracted_counts < 0] = 0
    if plot:
        plt.step(centers, subtracted_counts, label = 'subtracted hist')
        plt.legend()
        
    big_n = np.sum(subtracted_counts)
    norm_subtract_hist = subtracted_counts / big_n
    # weights = 1.0 / subtracted_counts / 
    mean_value = np.sum(centers * norm_subtract_hist)
    if plot:
        plt.title(f'Bias: {bias}V, mean_value: {mean_value}')
        plt.axvline(x = mean_value, color = 'green')  
    # mean_err = np.sum((centers/big_n) ** 2)(subtracted_counts) + (np.sum(subtracted_counts*centers)/(big_n)**2) ** 2 * (np.sum(subtracted_counts)) #overestimation
    a = np.sum(subtracted_counts * centers)
    mean_err = np.sqrt(np.sum( ((a - centers * big_n)/ big_n ** 2) ** 2 * (counts1 + counts2)))
    return (mean_value, mean_err)
#%% in this case the subtracted histogram is not useful since the LED wasn't bright enough. run cell to inspect it
#for run in runs:
#    get_subtract_hist_mean(run.bias,run.all_led_peak_data, run.all_dark_peak_data, plot = True)
#%% quickly get the approximate locations of the SPE peaks
cutoffs = []
centers_list=[]
centers_guesses = []
numpeaks = []
for run in runs:
    bins = int(round(np.sqrt(len(run.all_peak_data))))
    count, edges = np.histogram(run.all_peak_data, bins=bins)
    centers = (edges[:-1] + edges[1:])/2
    peaks, props = signal.find_peaks(count, prominence=22, distance=6) # prominence and distance are based on the particular properties of this data set and aren't always the same! read scipy documentation
    #print('peaks: ', peaks)
    high = min(len(peaks)-1,3)
    numpeaks.append(high)
    
    plt.figure()
    plt.hist(run.all_peak_data, bins=bins)
    #fig.savefig(f'bias = {run.bias}.png')
    
    fitrange = ((centers[peaks[high]] - centers[peaks[0]])/2)
    range_low =  centers[peaks[0]]- 0.25*fitrange
    range_high = centers[peaks[high]]+ 0.35*fitrange
    
    cutoffs.append((range_low, range_high))
    centers_list.append(centers[peaks[0]])
    peaks = peaks[0:]
    centers_guesses.append([centers[peak] for peak in peaks])
    
    #fig = plt.figure()
    #plt.hist(run.all_peak_data, bins=bins)
    #for peak in peaks:
    #    plt.scatter(centers[peak],count[peak])
    #plt.axvline(range_low, c='red')
    #plt.axvline(range_high, c='black')
    ## plt.xlim([0,1])
    #plt.yscale('log')
    #fig.savefig(f'bias = {run.bias}.png')
#%% testing cell - plot with peak fit. in this case the settings worked for all voltages
# set conditions, temp
#n=3
#T = 170
#con = 'LXe'
#info_spe = MeasurementInfo()
#info_spe.condition = con
#info_spe.date = runs[n].date
#info_spe.temperature = T
#info_spe.bias = biases[n]
#info_spe.baseline_numbins = 50
#info_spe.peaks_numbins = 200
#info_spe.data_type = 'h5'
#wp = WaveformProcessor(info_spe, run_info_self = runs[n], run_info_solicit = pedestal, baseline_correct = True, cutoff = cutoffs[n], centers = centers_guesses[n], peak_range = (1,numpeaks[n]+1), prom = prom)
#wp.process(do_spe = True, do_alpha = False)
#wp.plot_peak_histograms(log_scale = False)
#wp.plot_spe()
#%% plot all the fits and save them
plt.close('all')

savepath = 'plots/'
campaign_spe = []
for i in range(len(runs)):
    info_spe = MeasurementInfo()
    info_spe.condition = 'LXe'
    info_spe.date = runs[i].date
    info_spe.temperature = 169.6
    info_spe.bias = runs[i].bias
    info_spe.baseline_numbins = 50
    info_spe.peaks_numbins = 200
    info_spe.data_type = 'h5'
    wp_spe = WaveformProcessor(info_spe, run_info_self = runs[i], run_info_solicit = pedestal, baseline_correct = True,  cutoff = cutoffs[i], centers = centers_guesses[i], peak_range = (1,numpeaks[i]+1), prom = mprom * ( info_spe.bias - 32 ) + bprom)
    print('beginning wp_spe.process...')
    wp_spe.process(do_spe = True, do_alpha = False)
    wp_spe.plot_peak_histograms(log_scale = False, savefig=True, path = savepath+'gauss_'+str(info_spe.bias)+'.png')
    wp_spe.plot_spe(savefig=True, path = savepath+'spe_'+str(info_spe.bias)+'.png')
    #wp_spe.plot_both_histograms(log_scale = True, with_fit=True, with_baseline_fit=True, savefig=False, path = savepath+'baseline_'+str(info_spe.bias)+'.png')
    wp_spe.plot_peak_histograms(log_scale = False, savefig=True, path = savepath+'gauss_'+str(info_spe.bias)+'.svg')
    wp_spe.plot_spe(savefig=True, path = savepath+'spe_'+str(info_spe.bias)+'.svg')
    #wp_spe.plot_both_histograms(log_scale = True, with_fit=True, with_baseline_fit=True, savefig=False, path = savepath+'baseline_'+str(info_spe.bias)+'.svg')
    #plt.show()
    if wp_spe.get_spe()[1] != 0:
        campaign_spe.append(wp_spe)
    else:
        print(f"The following bias voltage has zero standard deviation in the SPE amplitude: {wp_spe.info.bias} V.\nThe SPE amplitude claims to be exactly {wp_spe.get_spe()[0]}")
    #plt.show()
#%% plot linear fit to the breakdown voltage
curr_campaign = campaign_spe

#wp = campaign_spe[1]
#spe_vals = wp.get_spe()[0]
#absolute_spe_vals = spe_vals / (invC_filter * 1.60217662e-7)
#spe_err = wp.get_spe()[1]

#print('absolute spe vals: ',absolute_spe_vals)
#print('absolute spe err: ',absolute_spe_vals * np.sqrt(
#                (spe_err * spe_err) / (spe_vals * spe_vals)
#                + (invC_err_filter * invC_err_filter) / (invC_filter * invC_filter)))
#
#print('cutoffs: ', cutoffs[1])
#print('center guesses: ', centers_guesses[1])
#print('numpeaks: ', numpeaks[1])
spe = SPE_data(curr_campaign, invC_filter, invC_err_filter, filtered = True, prominence = (mprom,bprom))
spe.plot_spe(in_ov = False, absolute = False,savefile=f'{n},{m}.svg')
#slop = spe.absolute_spe_res.params['slope'].value
#sloperr = spe.absolute_spe_res.params['slope'].stderr
#inter = spe.absolute_spe_res.params['intercept'].value
#intererr = spe.absolute_spe_res.params['intercept'].stderr
#redchi = spe.absolute_spe_res.redchi
#print(spe.spe_res.params["slope"].value)
#print(spe.spe_res.params["intercept"].value)
#print(spe.v_bd)

#plt.show()