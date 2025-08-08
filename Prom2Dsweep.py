# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from RunInfo import RunInfo
import matplotlib.pyplot as plt
from ProcessWaveforms_MultiGaussian import WaveformProcessor as WaveformProcessor
import os
import math
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from multiprocessing import Process,Queue
import time
from sweep_definitions import do_analysis,sendemail


if __name__ == '__main__':
    sendemail(subject = 'STARTED PROMINENCE SWEEP',contents ='',filename = 'prominence_sweep_emails.txt')
    time_start = time.time() #for bookkeeping

    #manager = Manager()
    Q = Queue()

    #1us, 10x gain, filter on - see Noah's report on OneDrive
    invC_filter = 0.01005
    invC_err_filter = 0.0005
    # 
    LOWERPROM = 0.0045 #Smallest prominence allowed
    UPPERPROM = 0.0065 #Largest prominence allowed
    PINTERCEPT = np.linspace(LOWERPROM,UPPERPROM,50) #creates a numpy array of prominence intercept values. These are the prominence values at bias = 32 V
    step = 0.00002 #step size for the slope of the prominence lines
    Z = np.full(shape = (len(PINTERCEPT),math.ceil((PINTERCEPT[-1]-PINTERCEPT[0])/(step*4))), fill_value=-1.0,dtype=np.float64) #setup for ult

    ult = dict({'Slope': Z.copy(), 'Slope error': Z.copy(), 'Intercept': Z.copy(), 'Intercept Error': Z.copy(), 'Chi-squared': Z.copy(), 'Breakdown': Z.copy(), 'Breakdown Error': Z.copy()}) #ult is the dictionary that stores all values once they are calculated

    path ='20250708_SPE_LED/' #path to the folder with the hdf5 run files
    date = '20250806'
    savepath = 'Sweep Lines/'+path+date+'/'
    title = date+' Prominence Sweep; Data = July 8, 2025 LED SPE' #also the filename for storing data
    upperlim_overrides = {32:.119}
    centers_overrides = {32:[.005, .009, .013],
                        32.5:[.0055, .01, .015, .021],
                        33:[0.006, 0.011, 0.017, 0.023],
                        33.5:[.0065,.013,.018,.024],
                        34:[0.007, 0.013, 0.020, 0.026],
                        34.5: [0.007,0.014, 0.021, 0.028],
                        35:[0.0075,0.015,0.022,0.030],
                        35.5:[0.008,0.016, 0.023,0.031],
                        36:[0.0085,0.017, 0.025, 0.033]
                        }

    #  baseline noise - 20V
    f = 'Run_1751987061.hdf5' #path to the baseline noise/"pedestal"
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
    # 
    # in general upper_limit is best set at scoperange - offset - 0.001; not always
    
    #upperlim = upperlim[:-1] #was causing problems so removed this. Usually, this is here to get rid of the value corresponding to the baseline

    for i in range(len(PINTERCEPT)): #loop through the intercepts of the prominence lines
        time_process_start = time.time()
        bprom = PINTERCEPT[i] #grab bias = 32V intercept for current prominence line
        PSLOPE = np.arange(0,(PINTERCEPT[-1]-bprom)/4,step) #create an array of slopes for this line. The slopes will never result in a prominence higher than UPPERPROM
        linedatas = [[bprom,i,PSLOPE[j],j] for j in range(len(PSLOPE))]
        processlist = [Process(target=do_analysis,args=(linedata,
                                                        biases,
                                                        scoperange,
                                                        offset,
                                                        pedestal,
                                                        path,
                                                        files,
                                                        [invC_filter,invC_err_filter],
                                                        Q,
                                                        savepath+f'{linedata[1]},{linedata[3]}/',
                                                        upperlim_overrides,
                                                        centers_overrides)) for linedata in linedatas]
        N = 5 #how many processes? How many runs to do at the same time?
        processlistlist = [processlist[i:i+N] for i in range(0,len(processlist),N)] #get a list of list of processes
        for l in processlistlist: #compute one list of analyses within processlistlist at a time until all processes are done
            [p.start() for p in l]
            [p.join() for p in l]
        for i in range(len(processlist)): #for every process scheduled, try to get an output
            out = Q.get() #get the output of the SPE analysis
            if out: #if there is one, record it.
                n = out[0]
                m = out[1]
                ult['Slope'][n][m] = out[2]
                ult['Slope error'][n][m] = out[3]
                ult['Intercept'][n][m] = out[4]
                ult['Intercept Error'][n][m] = out[5]
                ult['Chi-squared'][n][m] = out[6]
                ult['Breakdown'][n][m] = out[7]
                ult['Breakdown Error'][n][m] = out[8]
                #sendemail(subject = f'FINISHED PROMINENCE LINE: {n},{m}',contents = str(out),filename = 'prominence_sweep_emails.txt')
        #print(ult['Chi-squared'])
        time_process_end = time.time()
        time_process_total = time_process_end - time_process_start
        if len(processlist) != 0: # can't divide by zero
            rate_process = time_process_total/len(processlist)
        else:
            rate_process = time_process_total #need a number
        sendemail(subject='PROMINENCE SWEEP UPDATE',contents=f'Total time: {time_process_total/60:0.4} minutes. \n{rate_process:0.4} seconds per run. \n{60/rate_process:0.4} runs per minute. \n'+np.array2string(ult['Breakdown']),filename = 'prominence_sweep_emails.txt')
        
    #grab the stuff from ult
    slope = np.array(ult['Slope'])
    slope_err = np.array(ult['Slope error'])
    intercept = np.array(ult['Intercept'])
    intercept_err = np.array(ult['Intercept Error'])
    redchi = np.array(ult['Chi-squared'])
    breakdown = np.array(ult['Breakdown'])
    breakdown_err = np.array(ult['Breakdown Error'])
    logredchi = np.log10(redchi, out=np.full_like(redchi, fill_value=-10,dtype = float), where=(redchi>=0)) #for log scale stuff

    #setup for later stuff

    df_slope = pd.DataFrame( slope,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_slope.to_csv(title+'; slope.txt') #save the slope data
    df_slope_err = pd.DataFrame( slope_err,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_slope_err.to_csv(title+'; slope_err.txt') #save the slope error data
    df_intercept = pd.DataFrame( intercept,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_intercept.to_csv(title+'; intercept.txt') #save the itnercept data
    df_intercept_err = pd.DataFrame( intercept_err,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_intercept_err.to_csv(title+'; intercept_err.txt') #save the interecpt error data
    df_redchi= pd.DataFrame(  redchi,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_redchi.to_csv(title+'; redchi.txt') #save the reduced chi-squared data
    df_breakdown= pd.DataFrame(  breakdown,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_breakdown.to_csv(title+'; breakdown.txt') #save the breakdown data
    df_breakdown_err= pd.DataFrame(  breakdown_err,index = PINTERCEPT,columns = np.arange(0,(PINTERCEPT[-1]-PINTERCEPT[0])/4,step))
    df_breakdown_err.to_csv(title+'; breakdown_err.txt') #save the breakdown error data

    plt.close('all') #close all figures that might be open


    #plotting the chi-squares in a heatmap. 
    #yellow = large chi-square. black = small chi-square. Want close to 1.
    #chi-square = -1 is no data
    figdark,axdark = plt.subplots(1,1)

    dark = axdark.imshow(redchi, cmap = 'hot', interpolation = 'none', aspect='auto',extent=[0,(UPPERPROM-LOWERPROM)/4,LOWERPROM,UPPERPROM], origin='lower')
    axdark.set_xlabel("slope")
    axdark.set_ylabel('intercept')
    figdark.colorbar(dark)
    figdark.savefig(date+'_2dsweep_preliminary_darkzero.png')
    plt.close()


    #plotting the log of the chi-squares in a heatmap. 
    #yellow = large chi-square. black = small chi-square. Want close to 0.
    #log chi-square = -10 is no data
    figlogdark,axlogdark = plt.subplots(1,1)

    logdark = axlogdark.imshow(logredchi, cmap = 'hot', interpolation = 'none', aspect='auto',extent=[0,(UPPERPROM-LOWERPROM)/4,LOWERPROM,UPPERPROM], origin='lower')
    axlogdark.set_xlabel("slope")
    axlogdark.set_ylabel('intercept')
    figlogdark.colorbar(logdark)
    figlogdark.savefig(date+'_2dsweep_preliminary_darkzero_log.png')
    plt.close()



    #plotting the chi-squares in a heatmap. 
    #yellow = small chi-square. black = large chi-square. Want close to -1.
    #chi-square = 1 is no data
    figlit,axlit = plt.subplots(1,1)

    lit = axlit.imshow(-1*redchi, cmap = 'hot', interpolation = 'none', aspect='auto',extent=[0,(UPPERPROM-LOWERPROM)/4,LOWERPROM,UPPERPROM], origin='lower')
    axlit.set_xlabel("slope")
    axlit.set_ylabel('intercept')
    figlit.colorbar(lit)
    figlit.savefig(date+'_2dsweep_preliminary_lightzero.png')
    plt.close()



    #plotting the log of the chi-squares in a heatmap. 
    #yellow = small chi-square. black = large chi-square. Want close to 0.
    #chi-square = 10 is no data
    figloglit,axloglit = plt.subplots(1,1)

    loglit = axloglit.imshow(-1*logredchi, cmap = 'hot', interpolation = 'none', aspect='auto',extent=[0,(UPPERPROM-LOWERPROM)/4,LOWERPROM,UPPERPROM], origin='lower')
    axloglit.set_xlabel("slope")
    axloglit.set_ylabel('intercept')
    figloglit.colorbar(loglit)
    figloglit.savefig(date+'_2dsweep_preliminary_lightzero_log.png')
    plt.close()



    # use a different method of taking the logarithm for the fancier plots (give the no-data pixels an imaginary value)
    logredchi = np.log10(redchi, out=np.full(shape=redchi.shape, fill_value = -1j), where=(redchi>=0))


    #plotting the chi-squares in a heatmap. 
    #green = large chi-square. light blue = small chi-square. Want close to 1.
    #chi-square = -1 is no data
    figcut,axcut = plt.subplots(1,1,layout = 'constrained')

    normcut = colors.Normalize(vmin=-1,vmax=2,clip=True)
    cut = axcut.imshow(redchi, cmap = 'inferno',norm = normcut, interpolation = 'none', aspect='auto',extent=[0,(UPPERPROM-LOWERPROM)/4,LOWERPROM,UPPERPROM], origin='lower')
    axcut.set_xlabel("slope")
    axcut.set_ylabel('intercept')
    figcut.colorbar(cut)
    cutcolor = [cut.cmap(cut.norm(-1)),cut.cmap(cut.norm(2))]
    cutpatch = [mpatches.Patch(color = cutcolor[1], label = "reduced $\chi^2 \geq 2$"),mpatches.Patch(color = cutcolor[0], label = "No Data")]
    figcut.legend(handles=cutpatch,loc='outside right upper')
    #axcut.set_title()
    figcut.savefig(date+'_2dsweep_preliminary_cutoff_linear.png')
    plt.close()




    #plotting the log of the chi-squares in a heatmap. 
    #white = large chi-square. grey = small chi-square. black = no data
    figbinary,axbinary = plt.subplots(1,1,layout = 'constrained')

    greatermask = (logredchi.real >= 0.5)
    mediummask = (logredchi.real < 0.5) & (logredchi.imag == 0)
    lessermask = (logredchi.imag != 0)
    binarry = logredchi.copy()
    binarry[greatermask] = 1
    binarry[mediummask] = 0
    binarry[lessermask] = -1
    binarry = binarry.astype(dtype = int)
    zero_one = [-1,0,1]
    normbinary = colors.Normalize(vmin=-0.5,vmax=1,clip=True)
    binary = axbinary.imshow(binarry, cmap = 'gray',norm = normbinary, interpolation = 'none', aspect='auto',extent=[0,(UPPERPROM-LOWERPROM)/4,LOWERPROM,UPPERPROM], origin='lower')
    axbinary.set_xlabel("slope")
    axbinary.set_ylabel('intercept')
    blackwhite = [binary.cmap(binary.norm(i)) for i in zero_one]
    binarypatches = [mpatches.Patch(color=blackwhite[2], label="$log_{10}(\chi^2) \geq 0.5$"), mpatches.Patch(color=blackwhite[1], label="$log_{10}(\chi^2) < 0.5$"), mpatches.Patch(color=blackwhite[0], label="No Data")]
    figbinary.legend(handles=binarypatches,loc='outside right upper')
    #axbinary.set_title()
    figbinary.savefig(date+'_2dsweep_preliminary_binary_log.png')
    plt.close()

    number_of_runs = np.count_nonzero(ult['Chi-squared']+1)
    time_end = time.time()
    time_total = (time_end-time_start)
    rate = time_total/number_of_runs
    content = ('The sweep named "'+title+'" is complete.'
               +f'\nTotal time: {time_total/3600:0.4} hours'
               +f'\n{rate:0.4} seconds per run. {60/rate:0.4} runs per minute.'
               +'\n'+np.array2string(ult['Chi-squared']))
    sendemail(subject='SWEEP COMPLETE.',contents=content,filename="prominence_sweep_emails.txt")
    print("ALL DONE")