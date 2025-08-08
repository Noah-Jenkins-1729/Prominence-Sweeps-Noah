# Prominence-Sweeps-Noah
Prominence relationship analysis via parameter sweep for the SPE Analysis within PocarLab.

# Overview of Files
AnalyzePDE: This holds the SPEData class for making and fitting the SPE Amplitude plots.

MeasurmentInfo: This holds the MeasurementInfo class that is not used for anything. It is vestigial from past ways.

nexo_new.mplstyle: The matplotlib style configuration for standardized plot styles

ProcessWaveforms_hacked: The same as ProcessWaveforms_MultiGaussian, but fits an exponential with to the baseline of the fingerplots instead of a line.

ProcessWaveforms_MultiGaussian: Holds the WaveformProcessor class that creates and manages the fingerplot of a specific run with a linear fit for the baseline.

Prom2Dsweep: this is the main prominence sweep code. More on this later.

prominence_sweep_emails: a text file that the sendemail function in sweep_definitions

RunAnalysis: the standard run analysis code for doing an SPE analysis. almost nothing is automated or hidden within functions.

RunAnalysis_parasite: effectively the same as RunAnalysis, but it instead uses the same do_analysis function from sweep_definitions that Prom2Dsweep does

RunInfo: Holds the RunInfo class that processes the hdf5 file for a specific run. It extracts the peak heights with the given prominence value.

RunInfo_hacked: the same as RunInfo, except it grabs peak prominences instead of peak heights.

sweep_definitions: holds functions used for Prom2Dsweep and RunAnalysis_parasite for ease of changing parameters.

# Prom2Dsweep.py
To run a sweep, do the following:
1. Decide the range of intercepts you wish to take and input them into the UPPERPROM and LOWERPROM variables.
2. Decide how many rows you want on the heatmap and change the linspace for PINTERCEPT accordingly
3. decide the step size for the slope columns. If you want a number of columns, consult the definition of the Z ndarray, which has shape (#rows,#cols), and calculate the step size needed.
4. input the correct calibration constants for the invC_filter and invC_err_filter. These can be looked up in the Electronics Calibration report.
5. specify the path to the folder with the hdf5 files using the variable path. additionally make sure that your hdf5 files exist there and are uncorrupted/openable
6. set the current date in the date variable
7. decide if you want to save fingerplots and SPE Amplitude plots for each pixel. If no, go to line 94 (where savepath is called again) and input None. If yes, set the savepath variable to the path of the folder where you want to store the plots. Additionally, check the sweep_definitions file for where it checks the savefile variable and make sure the code there is either commented our or uncommented out (depending on your preference).
8. set the title variable to what you want to name the sweep. Typically, you'll want to put what dataset you're using in the name. It will automatically append the date to the front of the title.
9. set upperlim_overrides dictionary to whatever the optimal upperlims are for each bias voltage. You do not need to set anything, and you don't need to set all of them if you only want to set one. Upperlim sets an upperbound for how large peak heights can get.
10. repeat step six, but for the centers_overrides. Instead of floats, put lists for where the WaveformProcessors should guess where the centers of the gaussians are on the fingerplot.
11. ensure the baseline data is set correctly in the f variable. it should be labeled with a comment.
12. Decide how many cores you wish to be using at any given moment. It is usually a good idea not to exceed the number of cores your machine has. Input that in the N variable on line 97, just after where the processlist is defined but before the processlistlist definition.
13. if you want to get email updates as the sweep progresses, find all the sendemail function definition in sweep_definitions and give it your email address and an app password. You can also adjust what content you get in the emails by finding the sendemail function calls in the other sweep_definitions functions and in Prom2Dsweep and adjusting accordingly, or even commenting out the ones you don't want.
14. Prom2Dsweep will make some plots automatically. If you wish to adjust what it makes, or if you don't want it to make anything, you can find that code below the dataframe saving block, after the plt.close('all') function call. It is perfectly acceptable to remove them, as all the data used to make them is saved into csv files through the pandas.DataFrame class's to_csv method.
15. you should now be ready to run the sweep. run Prom2Dsweep.py however you wish, and wait for it to finish
