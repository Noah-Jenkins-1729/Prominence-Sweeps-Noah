import numpy as np
from MeasurementInfo import MeasurementInfo
#from RunInfo import RunInfo
from RunInfo_hacked import RunInfo_hacked as RunInfo
from scipy import signal
from AnalyzePDE import SPE_data
#from ProcessWaveforms_MultiGaussian import WaveformProcessor as WaveformProcessor
from ProcessWaveforms_hacked import WaveformProcessor_hacked as WaveformProcessor
import smtplib
from email.message import EmailMessage
from pathlib import Path
from multiprocessing import Queue

#the following function will email me with the contents string
def sendemail(subject,contents,filename):
    with open(filename,'w') as f: #filename must already exist
        f.write(contents) #writes contents into filename

    with open(filename) as fp:
        msg = EmailMessage() #create EmailMessage object
        msg.set_content(fp.read()) #sets the contents of msg to be the contents of filename

    msg['Subject'] = subject #subject of the email
    msg['From'] = "e@mail.ext" #from me
    msg['To'] = "e@mail.ext" #to me

    gmail = smtplib.SMTP('smtp.gmail.com',587) #standard for gmail
    gmail.ehlo() #no idea
    gmail.starttls() #no idea
    gmail.login("e@mail.ext","APP PASSWORD") #first entry is username, second is password. The second entry is an app password (2-Factor authentication prevents using standard password)
    gmail.send_message(msg) #send the email


def calculate(data):
    file = data[0]
    upperlim = data[1]
    prom = data[2]
    path = data[3]
    files = data[4]
    out = RunInfo([path+files[file]], specifyAcquisition = False, do_filter = True, upper_limit = upperlim, poly_correct=True, baseline_correct = True, prominence = prom, is_led = True)
    return out

def histograms(data):
    #try:
    run = data[0]
    centers_override = data[1]

    bins = int(round(np.sqrt(len(run.all_peak_data))))
    count, edges = np.histogram(run.all_peak_data, bins=bins)
    centers = (edges[:-1] + edges[1:])/2
    peaks = None
    if centers_override:
        for i in centers_override.keys():
            if i == run.bias:
                peaks = [np.argmin(np.abs(np.array(centers)-peak)) for peak in centers_override[i]]
    if not peaks:
        peaks, props = signal.find_peaks(count, prominence=22, distance=6) # prominence and distance are based on the particular properties of this data set and aren't always the same! read scipy documentation
    high = min(len(peaks)-1,3)
    numpeak = high

    #plt.figure()
    #plt.hist(run.all_peak_data, bins=bins)
    
    fitrange = ((centers[peaks[high]] - centers[peaks[0]])/2)
    range_low =  centers[peaks[0]]- 0.25*fitrange
    range_high = centers[peaks[high]]+ 0.35*fitrange
    
    cutoff = (range_low, range_high)
    center_list = (centers[peaks[0]])
    peaks = peaks[0:]
    center_guesses = [centers[peak] for peak in peaks]
    #plt.close()
    return [cutoff,center_list,center_guesses,numpeak]
    #except:
    #    sendemail(subject = f'ALERT: ERROR DURING HISTOGRAM CALCULATION', contents = f'Histogram failed on bias {run.bias} V', filename="prominence_sweep_emails.txt")
    #    return None

def make_campaign(data: tuple[RunInfo,int,float,tuple[float,float],tuple,list[float],int,RunInfo,str]
                  ):
    run = data[0]
    i = data[1]
    prom = data[2]
    cutoff = data[3]
    center_list = data[4]
    center_guess = data[5]
    numpeak = data[6]
    pedestal = data[7]
    savepath = data[8]
    info_spe = MeasurementInfo() #initializing lines...
    info_spe.condition = 'LXe'
    info_spe.date = run.date
    info_spe.temperature = 169.6
    info_spe.bias = run.bias
    info_spe.baseline_numbins = 50
    info_spe.peaks_numbins = 200
    info_spe.data_type = 'h5'

    #try: #try is for debugging purposes and to prevent any errors from hindering the sweep
    wp_spe = WaveformProcessor(info_spe, run_info_self = run, run_info_solicit = pedestal, baseline_correct = True,  cutoff = cutoff, centers = center_guess, peak_range = (1,numpeak+1), prom = prom)
    wp_spe.process(do_spe = True, do_alpha = False) #calcualte the SPE amplitudes for each bias voltage
    if savepath:
        #Path(savepath).mkdir(parents=True,exist_ok=True)
        #wp_spe.plot_peak_histograms(log_scale = False, savefig=True, path = savepath+'gauss_'+str(info_spe.bias)+'.png')
        #wp_spe.plot_spe(savefig=True, path = savepath+'spe_'+str(info_spe.bias)+'.png')
        pass
    #except:
    #    sendemail(subject = f'ALERT: ERROR DURING WAVEFORMPROCESSOR CALCULATION', contents = f'WaveformProcessor failed on bias {run.bias} V with prominence {prom} V', filename="prominence_sweep_emails.txt")
    #    return None
    if (wp_spe.get_spe()[1] != 0) and (not (np.array(wp_spe.peak_stds)==None).any()):
        return wp_spe #if the data is correct, add the WaveformProcessor to the campaign; it is ready to be fed into SPEData
    else:
        sendemail(subject = f'ALERT: ERROR DURING WAVEFORMPROCESSOR CALCULATION', contents = f'WaveformProcessor failed on bias {run.bias} V with prominence {prom} V', filename="prominence_sweep_emails.txt")
        return None


def do_analysis(linedata: tuple[float,int,float,int],
                bias_list: list[float],
                scoperange_list: list[float],
                offset_list: list[float],
                pedestal: RunInfo,
                path: str,
                files_list: list[str],
                invC: tuple[float,float],
                Q: Queue,
                savepath: str = None,
                upperlim_override: dict = None,
                centers_override: dict = None
                ):
    bprom = linedata[0]
    n = linedata[1]
    mprom = linedata[2]
    m = linedata[3]
    biases = np.array(bias_list)
    scoperange = scoperange_list
    offset = offset_list
    pedestal = pedestal
    path = path
    files = files_list
    invC_filter = invC[0]
    invC_err_filter = invC[1]

    upperlim = np.array(scoperange)-np.array(offset)-0.001
    if upperlim_override:
        for i in upperlim_override.keys():
            upperlim = [upperlim_override[i] if v==i else u for v,u in zip(biases,upperlim)]
    
    #print('\n')
    #print(biases)
    #print(upperlim)
    #input('Press enter to continue...')

    print(f'\nCURRENT PROMINENCE LINE: {mprom} * (bias - 32) + {bprom}')
    f_upperlim_prom_path_file = [[f,upperlim[f],mprom * ( biases[f] - 32 ) + bprom,path,files] for f in range(len(biases)) if biases[f]>25]
    res = map(calculate,f_upperlim_prom_path_file)
    runs = list(res)
    runs = [run for run in runs if run != None]
    runs = [run for run in runs if run.bias > 25]
    res = map(histograms,[[run,centers_override] for run in runs])
    hists = list(res)
    hists = [hist for hist in hists if hist != None]
    cutoffs = [out[0] for out in hists]
    centers_list = [out[1] for out in hists]
    centers_guesses = [out[2] for out in hists]
    numpeaks = [out[3] for out in hists]
    run_f_prom_cutoff_center_guess_numpeak_pedestal_save_cover = [[runs[i],f_upperlim_prom_path_file[i][0],f_upperlim_prom_path_file[i][2],cutoffs[i],centers_list[i],centers_guesses[i],numpeaks[i],pedestal,savepath,centers_override] for i in range(len(runs))]
    res = map(make_campaign,run_f_prom_cutoff_center_guess_numpeak_pedestal_save_cover)
    curr_campaign = list(res)
    curr_campaign = [l for l in curr_campaign if l != None]

    try: #the SPEData object will calculate the SPE amplitude vs bias voltage relationship
        spe = SPE_data(curr_campaign, invC_filter, invC_err_filter, filtered = True, prominence = [mprom,bprom])
        if savepath:
            #Path(savepath).mkdir(parents=True,exist_ok=True)
            #spe.plot_spe(savefile = savepath+f'{n},{m}.png')
            pass
        bit = [n,m,
            spe.spe_res.params["slope"].value,
            spe.spe_res.params["slope"].stderr,
            spe.spe_res.params["intercept"].value,
            spe.spe_res.params["intercept"].stderr,
            spe.spe_res.redchi,
            spe.v_bd,
            spe.v_bd_err
            ]
        Q.put(bit)
        print(f'FINISHED PROMINENCE LINE: {mprom} * (bias - 32) + {bprom}\n')
        #sendemail(subject = f'FINISHED PROMINENCE LINE: {mprom} * (bias - 32) + {bprom}',contents = str(bit),filename = 'prominence_sweep_emails.txt')
    except: #if there is an error, skip this prominence line, record it, and send an email
        sendemail(subject = f'ALERT: PROMINENCE LINE FAILED',contents = f'SPEData failed with prominence line {mprom} * (bias - 32) + {bprom}', filename="prominence_sweep_emails.txt")
        Q.put(False)
