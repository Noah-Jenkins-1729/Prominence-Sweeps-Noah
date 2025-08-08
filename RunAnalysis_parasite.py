from sweep_definitions import do_analysis
import numpy as np
import pandas as pd
from multiprocessing import Queue
import os
from RunInfo import RunInfo

df_redchi = pd.read_csv('20250715 Prominence Sweep; Data = October 17, 2024 LED SPE; redchi.txt',index_col=0)
PINTERCEPT = np.array(df_redchi.index.values)
PSLOPE = np.array(df_redchi.columns.astype('float').values)

n = 0
m = 0
bprom = PINTERCEPT[n] #intercept at bias = 32
mprom = PSLOPE[m]  #slope of prominence line


path ='20241017_SPE_LED_405_LXe/' #path to the folder with the hdf5 run files
savefile = f'Sweep Lines/20241017_SPE_LED_405_LXe/{n},{m}/'

#  baseline noise - 20V
f = 'Run_1729184064.hdf5' #path to the baseline noise/"pedestal"
pedestal = RunInfo([path+f],  specifyAcquisition = False, do_filter = True, is_solicit = True, upper_limit = 0.008, baseline_correct = True,  poly_correct=False)

# loop quickly gets the file names
# separate files for each bias voltage -> specifyAcquisition = False
files = os.listdir(path)
# Filtering only the files
files = [f for f in files if os.path.isfile(path+'/'+f)]
#print(files)
#input("Press enter to continue...") #this and the previous line are useful for checking what files are being grabbed

#  read minimal waveforms as a shortcut to print all the metadata
runs = []
for file in range(len(files)):
    run_spe = RunInfo([path+files[file]], specifyAcquisition = False, do_filter = True, upper_limit = 5, poly_correct=True, baseline_correct = True, prominence = 0.008, is_led = True, num_waveforms = 1)
    runs.append(run_spe)
biases = list([run.bias for run in runs if run.bias > 25]) #grab biases list
#print(biases)
scoperange = list([run.yrange for run in runs if run.bias > 25])
#print(scoperange)
offset = list([run.offset for run in runs if run.bias > 25])
#print(offset)
files = list([files[f] for f in range(len(files)) if runs[f].bias > 25])

invC_filter = 0.0115
invC_err_filter = 0.000098
Q = Queue()

linedata = [bprom,n,mprom,m]
do_analysis(linedata,biases,scoperange,offset,pedestal,path,files,[invC_filter,invC_err_filter],Q,savefile)

out = Q.get()

print(out[7])
print(out[8])