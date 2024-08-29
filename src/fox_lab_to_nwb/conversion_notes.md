
# General comments
Point-person: Amy Streets (axs2909@case.edu) and Kris Lea (kxl786@case.edu)

https://sites.google.com/site/cwrufoxlab/people

# DeepLabCut output files
config files missing but Video and DLC output files are available

# Static TIFF files for confocal imaging
Not yet available

# Synchronization signal with Spike2
Not yet available

# Intracellular electrophysiology data including the voltage trace and stimulus trace

Not yet available

# Files on data shared 

### Tshx18D07_240124_115923_f3_r1.fly2
This is supposed to be a struct with DAQ data and the wingbeat analysis. Can't open it with matlab online unless extension is changed to .mat

Channel names:
- CamSync
- CamTrigger
- OptoTrigger
- LWingBeatAmp
- RWingBeatAmp
- WingBeatFreq
- LHutchen  # What does Hutchen mean? Could be halteres?
- RHutchen
- PTrigger

These can't be extracted with pyton from matlab because they are strings. Are they always the same order?

## Cameras 
### Fastec

Files:
* TOPCAM_000000.avi
* SIDE_000000.avi

Format is avi, metadata is in xml

https://www.fastecimaging.com/

TopCam has a DLC analysis for the following body parts:

bodyparts: L_ant, R_ant, L_base, R_base,

Sidecam does not have DLC analysis

### Phantom
Files:
* XZ_1_186.mp4

Format is mp4, metadata is in xml

This camera has an associated DLC analysis for the following body parts:

bodyparts: haltere

### Photron

Files:
* Tshx18D07_240124_115923_f3_r1_down20.mp4 [This seems to be the same as topcam]


Not sure about this. The `.mii` metadata is not available and there is this extra file:


It also does not have a corresponding DLC analysis


## `Tshx18D07_f3_r1_trial.Tsh`
This should be trial data.
Unable to read in matlab online unless extension is changed to .mat

## Tshx18D07_240124_115923_f3_r1_down20_PROC.mat
This has the analysis for the wingbeat redone.


## Tshx18D07_ant_top_f3_r1_500ms_dvProject.mat
Not described