# Ephys data

## Data organization

```bash
Command: 	 tree -L 2
.
├── data
│   ├── 2024_10_24_E2.mat  # This contains the trialized data
│   ├── DeepLabCut
│   ├── seal tests
│   └── videos
└── experiment scripts
    ├── ComputeCellStats.m
    ├── RunExperiment_OMN.m
    └── SaveSingleSealTest.m

6 directories, 4 files

```

The top level file looks like a trial table, then the videos contain videos for at least three angles in the shared that can be mapped to the trial table.


## The general metadata (top level file)

This is in a matlab file that is dataframe-like in its structure but some fields are 
arrays. To extract the data of the first line use the following code:

```python
from pymatreader import read_mat
matlab_data = read_mat(file_path)["data"]


non_array_keys = []
array_keys = []

for key, value in matlab_data.items():
    first_value = value[0]
    if isinstance(first_value, np.ndarray):
        array_keys.append(key)
        print(f"{key=} {first_value.shape=}")
    elif isinstance(first_value, list):
        print(f"{key=}, len(value)={len(value)}")
        array_keys.append(key)
    else:
        non_array_keys.append(key)
        print(key, first_value)
```

This should output something like this:


```
date 2024_10_24
expNumber 2
trial 1
condition antennae free
age 3-6dpe
genotype 24C06
samplerate 10000
fps 100
scaleCurrent 100
scaleVoltage 100
nframes 1000
key='Vm' first_value.shape=(100000,)
key='I' first_value.shape=(100000,)
key='filteredVm' first_value.shape=(100000,)
key='puffer' first_value.shape=(100000,)
key='tachometer' first_value.shape=(100000,)
```

What are the units of Vm, I, filteredVm, puffer, tachometer?


The non-list entries of this file can be thought as a data frame and you can get it like this:

```python
import pandas as pd 

non_array_data = {key: matlab_data[key] for key in non_array_keys}
df = pd.DataFrame(non_array_data)
df.head(n=15)
```
Which should output something like this.


|    |       date |   expNumber |   trial | condition     | age    | genotype   |   samplerate |   fps |   scaleCurrent |   scaleVoltage |   nframes |
|---:|-----------:|------------:|--------:|:--------------|:-------|:-----------|-------------:|------:|---------------:|---------------:|----------:|
|  0 | 2024_10_24 |           2 |       1 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  1 | 2024_10_24 |           2 |       2 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  2 | 2024_10_24 |           2 |       3 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  3 | 2024_10_24 |           2 |       4 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  4 | 2024_10_24 |           2 |       5 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  5 | 2024_10_24 |           2 |       6 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  6 | 2024_10_24 |           2 |       7 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  7 | 2024_10_24 |           2 |       8 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  8 | 2024_10_24 |           2 |       9 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
|  9 | 2024_10_24 |           2 |      10 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
| 10 | 2024_10_24 |           2 |      11 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
| 11 | 2024_10_24 |           2 |      12 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
| 12 | 2024_10_24 |           2 |      13 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
| 13 | 2024_10_24 |           2 |      14 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |
| 14 | 2024_10_24 |           2 |      15 | antennae free | 3-6dpe | 24C06      |        10000 |   100 |            100 |            100 |      1000 |




## Videos

They look like this:

```
2024_10_24_E2_Video_lateral_flyLeft1.avi
2024_10_24_E2_Video_lateral_flyLeft2.avi
2024_10_24_E2_Video_lateral_flyLeft3.avi
2024_10_24_E2_Video_lateral_flyLeft4.avi
2024_10_24_E2_Video_lateral_flyLeft5.avi
2024_10_24_E2_Video_lateral_flyLeft6.avi
2024_10_24_E2_Video_lateral_flyLeft7.avi
2024_10_24_E2_Video_lateral_flyLeft8.avi
```

And there are more videos for fly lateral right as well

```
2024_10_24_E2_Video_lateral_flyRight_1.avi
2024_10_24_E2_Video_lateral_flyRight_2.avi
2024_10_24_E2_Video_lateral_flyRight_3.avi
2024_10_24_E2_Video_lateral_flyRight_4.avi
```

And video lateral ventral:

```
2024_10_24_E2_Video_lateral_ventral_1.avi
2024_10_24_E2_Video_lateral_ventral_2.avi
2024_10_24_E2_Video_lateral_ventral_3.avi
2024_10_24_E2_Video_lateral_ventral_4.avi
2024_10_24_E2_Video_lateral_ventral_5.avi
```

The pattern is:
{date}_{expNumber}_{view}_{trial}.avi

And they can be mapped to the trials using the `trial` column in the top level file matlab file.


## Deep Lab Cut Data

It looks like this:

```
2024_10_24_E2_Video_lateral_flyLeft1DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000.h5
2024_10_24_E2_Video_lateral_flyLeft1DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000_meta.pickle
2024_10_24_E2_Video_lateral_flyLeft2DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000.h5
2024_10_24_E2_Video_lateral_flyLeft2DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000_meta.pickle
```

And it seems that it only contains pose estimation for the Lateral view.

## Seal test

```
2024_10_24_SealTest_14.mat
2024_10_24_SealTest_15.mat
2024_10_24_SealTest_16.mat
```

There is only three files, it seems that the the last number is the trial but it is not clear. 

The content of one of these fields is simpler and looks like:

```python
{'date': '2024_10_24',
 'expNumber': 1,
 'trial': 1,
 'numSecOut': 2,
 'scaleCurrent': 100,
 'scaleVoltage': 100,
 'SAMPLERATE': 10000,
 'Vm': array([5.62302898, 5.62302898, 5.62302898, ..., 5.59015469, 5.5572804 ,
        5.65590327], shape=(40000,)),
 'I': array([6.93800059, 6.97087488, 7.5626121 , ..., 6.74075485, 6.93800059,
        6.44488624], shape=(40000,))
}
```

And here is the table for the non-array values

|    | date   |   expNumber |   trial | notes   |   SAMPLERATE |   TRIAL_TIME_LIGHT |   PRE_TRIAL_TIME |   POST_TRIAL_TIME |   numSecOut | genotype   |   fps |   nframes |
|---:|:-------|------------:|--------:|:--------|-------------:|-------------------:|-----------------:|------------------:|------------:|:-----------|------:|----------:|
|  0 |        |           6 |       1 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  1 |        |           6 |       2 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  2 |        |           6 |       3 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  3 |        |           6 |       4 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  4 |        |           6 |       5 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  5 |        |           6 |       6 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  6 |        |           6 |       7 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  7 |        |           6 |       8 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  8 |        |           6 |       9 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
|  9 |        |           6 |      10 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
| 10 |        |           6 |      11 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
| 11 |        |           6 |      12 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
| 12 |        |           6 |      13 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
| 13 |        |           6 |      14 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |
| 14 |        |           6 |      15 | fly01   |        10000 |                  4 |                2 |                 4 |          10 | 24C06      |    60 |       600 |



## Patch clamp data questions
* For the patch clamp, why does the data says 1,000 frames but the length of the data is `100,000`
* samplerate: 10000, fps: 100 is the fps for the video?
* is dep days post eclosing (e.g. 6-8dpe) or days post patching?
* What is scale current?
* What is scale voltage?
* What does the tachometer measures or represents? It is usually used to measure the speed of a rotating object, such as a motor or other machine.
* What about the puffer?
* Description and units of `Vm`, `Im`, `tachometer`, `puffer`, `filteredVm`
* What device was used to record the data?
* Are all the experiments current clamp?
* Any information about the electrodes?


## Optogenetic Stimuli data
* `SAMPLERATE` (10_000) vs `fps` (60) vs nframes (6000)
* What are `TRIAL_TIME_LIGHT`, `PRE_TRIAL_TIME` and `POST_TRIAL_TIME`?
* A description and units of `stimTiming` and `micLeft`
* What are `numSecOut`?

It seems that `stimTiming` the optogenetic?
It seems that `notes` is the subject name


## Video

Both videos have:
Number of frames: 600
Frame rate: 60.0
Duration: 10.0

The filenames seem to be:
- `2023_11_30_E6_Video_lateral_10.avi`
- `{date}_{expNumber}_{view}_{trial}.avi`

# Behavioral Data (Mill's second batch)

## Data organization

```bash
pwd
/home/heberto/cohen_project/Sample data/Suver Lab/behavioralDataExamples
 tree -shL 2
[4.0K]  .
├── [ 80M]  2024_05_29_E3.mat  # Speedy Bars
├── [140M]  2024_07_08_E3.mat  # Windy Steps
├── [ 94M]  2024_10_28_E4.mat # Coco
├── [5.9K]  Coco.m
├── [161K]  readme.pdf
├── [4.5K]  speedyBars.m
└── [6.6K]  windySteps.m
```



## WindySteps Experiment
**File:** 2024_07_08_E3.mat

### Experimental Design
In the experiments shown in Figure 1-3 and Figure S1, each fly experienced 10 airflow trials either increasing from 0 cm/s to 300 cm/s (5 trials) or decreasing from 300 cm/s to 0 cm/s (5 trials). Airflow was presented in 50 cm/s steps over 42 seconds, with 6 seconds of steady state airflow for each windspeed. The order of increasing and decreasing trials was pseudo-randomized for each flight. For each windspeed, averages were taken over the last 3 seconds as the mass flow controller did not adjust instantaneously. 6 seconds of pre-trial at the initial windspeed (with the solenoid valve off, to prevent airflow) was fed to the mass flow controller to allow for adjustment and settlement at the correct speed before the trial began. Additionally, at 0 cm/s windspeed the solenoid valve was shut off to prevent airflow.

### Data Fields
- **date**: date of experiment
- **expnumber**: experiment number on above date
- **condition**: fly condition being run (refers to either genotype or if wild type, typically visual condition)
- **min_age**: minimum age of flies being run in this experiment
- **max_age**: maximum age of flies being run in this experiment
- **samplerate**: sampling rate (in Hz)
- **adjust_time**: adjustment time at each windspeed (i.e. time not used in average, explained in paragraph above)
- **record_time**: time at each windspeed used in averages later on (total time at each windspeed is adjust_time + record_time)
- **fps**: framerate of video recording
- **nframes**: total number of frames recorded
- **trial**: trial within current experiment
- **stimtype**: 1 for ascending/increasing stimuli, 2 for descending/decreasing stimuli
- **pufferSignal**: vector spanning the length of experiment (recorded at samplerate) of voltage reading from puffer. When a puff is applied, voltage spikes
- **tachometerSignal**: vector spanning the length of experiment (recorded at samplerate) of tachometer recording wingbeat data from fly. Recorded in volts
- **tachometerSignal_smoothed**: post-processed vector of tachometerSignal smoothed using a Butterworth bandpass filter (150-250 Hz) followed by a second-order Savitzky-Golay filter. Same length as tachometerSignal

### Additional Technical Details
From analysis of the `windySteps.m` MATLAB script:

#### Hardware Configuration
- Uses a National Instruments DAQ system (referenced as "Dev2")
- **Input Channels**:
  - Tachometer signal: analog input channel 0 (ai0)
  - Puffer signal: analog input channel 5 (ai5) with SingleEnded terminal configuration
- **Output Channels**:
  - Mass Flow Controller (MFC) signal: analog output channel 0 (ao0)
  - Optogenetic stimulation (when used): analog output channel 1 (ao1)
  - Solenoid valve control: analog output channel 2 (ao2)
  - Camera trigger: analog output channel 3 (ao3)

#### Wind Speed Calculation
- Wind tube radius: 0.18796 cm
- Wind speeds: 0, 50, 100, 150, 200, 250, 300 cm/s
- Wind speeds are converted from cm/s to voltage values for the MFC using the formula:
  ```
  tube_area = π × r² (in cm²)
  cm_flow_rate_per_sec = windspeeds × tube_area
  cm_flow_rate_per_min = cm_flow_rate_per_sec × 60
  liters_per_min = cm_flow_rate_per_min / 1000
  mfc_values = liters_per_min × 5/2
  ```

#### Signal Processing
- Tachometer signal is filtered using a 2nd-order Butterworth bandpass filter (150-250 Hz)
- Filter code: `[b,a] = butter(2,[150 250]/(Fs/2))`
- Filtered signal is further smoothed using Savitzky-Golay filtering

#### Video Recording
- Camera: Dorsal view (labeled as "Coronal (x) Camera")
- Resolution: 640x480 pixels, Y8 format
- Frame rate: 60 fps
- Manual shutter mode with shutter value of 800
- Videos saved as Motion JPEG AVI files
- Camera is hardware-triggered via the DAQ system

#### Experimental Conditions
- The script supports various experimental conditions including:
  - Normal flies
  - Silenced flies (with optogenetic manipulation)
  - Dark conditions
  - Various genetic lines (18D07, 74C10, ChrimCS, silencedCS, silencedCS_glued)

### How to read it:

```python
from pymatreader import read_mat
from pathlib import Path

folder_path = Path("/home/heberto/cohen_project/Sample data/Suver Lab/behavioralDataExamples")
file_path =  folder_path / "2024_07_08_E3.mat"
assert file_path.is_file()

mat_data = read_mat(file_path)["data"]
number_of_trials = len(mat_data)
print(f"number_of_trials={number_of_trials}")

non_array_keys = []
array_keys = []

for key, value in mat_data.items():
    first_value = value[0]
    if isinstance(first_value, np.ndarray):
        array_keys.append(key)
        print(f"{key=} {first_value.shape=}")
    elif isinstance(first_value, list):
        print(f"{key=}, len(value)={len(value)}")
        array_keys.append(key)
    else:
        non_array_keys.append(key)
        print(key, first_value)

number_of_trials=15
date None
expnumber 3
condition dark
min_age 3
max_age 5
samplerate 20000
adjust_time 3
record_time 3
fps 60
nframes 2520
trial 1
stimType 1
key='pufferSignal' first_value.shape=(960000,)
key='tachometerSignal' first_value.shape=(960000,)
key='tachometerSignal_smoothed' first_value.shape=(960000,)
```

## SpeedyBars Experiment
**File:** 2024_05_29_E3.mat

### Experimental Design
In the experiments shown in Figure 4, each fly experienced 60 optic flow trials ranging from 0–50 cm/s in 5 cm/s increments along with 100 cm/s (5 trials at each optic flow speed, pseudo-randomized order). Each trial consisted of an 8 second period of constant optic flow. For each optic flow speed, averages were taken over the last 6 seconds of each trial. The interlude between trial periods consisted of a non-moving static grating, equivalent to the 0 cm/s optic flow condition. The solenoid valve for airflow was shut off for the entirety of these experiments.

### Data Fields
- **date**: date of experiment
- **expnumber**: experiment number on above date
- **condition**: fly condition being run (refers to either genotype or if wild type, typically visual condition)
- **min_age**: minimum age of flies being run in this experiment
- **max_age**: maximum age of flies being run in this experiment
- **samplerate**: sampling rate (in Hz)
- **adjust_time**: adjustment time at each optic speed (i.e. time not used in average, explained in paragraph above)
- **record_time**: time at each windspeed used in averages later on (total time at each windspeed is adjust_time + record_time)
- **fps**: framerate of video recording
- **nframes**: total number of frames recorded
- **trial**: trial within current experiment
- **stimulus**: current optic flow speed being presented
- **pufferSignal**: vector spanning the length of experiment (recorded at samplerate) of voltage reading from puffer. When a puff is applied, voltage spikes
- **tachometerSignal**: vector spanning the length of experiment (recorded at samplerate) of tachometer recording wingbeat data from fly. Recorded in volts
- **photodiodeSignal**: vector spanning the length of experiment (recorded at samplerate) of a voltage reading from an amplified photodiode pointed at the visual stimuli. Signal increases when white bar is present relative to when black bar is present

### Additional Technical Details
From analysis of the `speedyBars.m` MATLAB script:

#### Hardware Configuration
- Uses a National Instruments DAQ system (referenced as "Dev2")
- **Input Channels**:
  - Tachometer signal: analog input channel 0 (ai0)
  - Photodiode signal: analog input channel 2 (ai2) with SingleEnded terminal configuration
  - Puffer signal: analog input channel 5 (ai5) with SingleEnded terminal configuration
- **Output Channels**:
  - Mass Flow Controller (MFC) signal: analog output channel 0 (ao0)
  - Solenoid valve control: analog output channel 2 (ao2)
  - Camera trigger: analog output channel 3 (ao3)

#### Visual Stimulation System
- Uses a TCP/IP server (port 5000) to communicate with a separate visual stimulus program (`videoClient.m`)
- Optic flow speeds: 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, and 100 cm/s
- 5 trials per stimulus speed, randomized order (60 trials total)
- Trial timing:
  - Adjustment time: 2 seconds
  - Recording time: 6 seconds
  - Total trial time: 8 seconds

#### Synchronization
- Uses a timer object for precise synchronization between visual stimuli and data acquisition
- 6-second delay between trial setup and start
- Sends start time to visual stimulus program via TCP/IP
- Waits for confirmation from visual stimulus program before proceeding

#### Video Recording
- Camera: Dorsal view (labeled as "Coronal (x) Camera")
- Resolution: 640x480 pixels, Y8 format
- Frame rate: 60 fps
- Manual shutter mode with shutter value of 800
- Videos saved as Motion JPEG AVI files
- Camera is hardware-triggered via the DAQ system

### How to read it:

```python
from pymatreader import read_mat
from pathlib import Path

folder_path = Path("/home/heberto/cohen_project/Sample data/Suver Lab/behavioralDataExamples")
file_path =  folder_path / "2024_05_29_E3.mat"
assert file_path.is_file()

mat_data = read_mat(file_path)["data"]
number_of_trials = len(mat_data)
print(f"number_of_trials={number_of_trials}")

non_array_keys = []
array_keys = []

for key, value in mat_data.items():
    first_value = value[0]
    if isinstance(first_value, np.ndarray):
        array_keys.append(key)
        print(f"{key=} {first_value.shape=}")
    elif isinstance(first_value, list):
        print(f"{key=}, len(value)={len(value)}")
        array_keys.append(key)
    else:
        non_array_keys.append(key)
        print(key, first_value)


number_of_trials=15
date None
expnumber 3
condition still
min_age 5
max_age 7
samplerate 20000
adjust_time 2
record_time 6
fps 60
nframes 480
trial 1
stimulus 35
key='pufferSignal' first_value.shape=(160000,)
key='tachometerSignal' first_value.shape=(160000,)
key='photodiodeSignal' first_value.shape=(160000,)
```

## Coco Experiment
**File:** 2024_10_28_E4.mat

### Experimental Design
In the experiments shown in Figure 5, each fly experienced 36 trials at 3 oscillatory frequencies (0.3, 1.3, 2.3 Hz, 12 trials at each frequency) with windspeeds ranging from 100–200 cm/s and optic flow speeds from 5–35 cm/s.

### Data Fields
- **date**: date of experiment
- **expnumber**: experiment number on above date
- **condition**: fly condition being run (refers to either genotype or if wild type, typically visual condition)
- **min_age**: minimum age of flies being run in this experiment
- **max_age**: maximum age of flies being run in this experiment
- **samplerate**: sampling rate (in Hz)
- **fps**: framerate of video recording
- **nframes**: total number of frames recorded
- **trialLength**: length of trial in seconds
- **block**: current block for current trial (experimental design consisted of each windspeed presented pseudo-randomly in a block, 12 blocks total were presented)
- **block_trial**: trial within the current block (either 1, 2, or 3)
- **trialnum**: trial within overall experiment
- **stimulus**: current oscillation frequency being presented
- **pufferSignal**: vector spanning the length of experiment (recorded at samplerate) of voltage reading from puffer. When a puff is applied, voltage spikes
- **tachometerSignal**: vector spanning the length of experiment (recorded at samplerate) of tachometer recording wingbeat data from fly. Recorded in volts
- **photodiodeSignal**: vector spanning the length of experiment (recorded at samplerate) of a voltage reading from an amplified photodiode pointed at the visual stimuli. Signal increases when white bar is present relative to when black bar is present

### Additional Technical Details
From analysis of the `Coco.m` MATLAB script:

#### Hardware Configuration
- Uses a National Instruments DAQ system (referenced as "Dev2")
- **Input Channels**:
  - Tachometer signal: analog input channel 0 (ai0)
  - Puffer signal: analog input channel 5 (ai5) with SingleEnded terminal configuration
  - Photodiode signal: analog input channel 6 (ai6) with SingleEnded terminal configuration
- **Output Channels**:
  - Mass Flow Controller (MFC) signal: analog output channel 0 (ao0)
  - Solenoid valve control: analog output channel 2 (ao2)
  - Camera trigger: analog output channel 3 (ao3)
- Uses external digital trigger (PFI0) with rising edge condition

#### Experimental Conditions
- Four possible conditions: 'wind', 'visual', 'both', or 'none'
- Condition determines whether wind oscillation, visual oscillation, both, or neither is applied

#### Oscillation Parameters
- Three oscillation frequencies: 0.3, 1.3, and 2.3 Hz
- Corresponding TCP frequencies for visual system: 3, 13, and 23
- Amplitude voltages: 0.835, 1.335, and 2.09 V
- Base voltage (150 cm/s wind speed): 2.4972571 V
- Wind oscillation formula: `((amplitude * sin(2*pi*t*frequency))+baseV)`

#### Trial Structure
- Total trial time: 19 seconds
- Video recording time: 16 seconds
- Baseline time: 6 seconds (constant wind speed)
- Oscillation time: 12 seconds
- End time: 1 second

#### Experimental Design
- 12 blocks with 3 trials per block (36 trials total)
- Each block contains all three frequencies in random order
- Block and trial-within-block structure is explicitly tracked

#### Visual Stimulation System
- Uses a TCP/IP server (port 5000) to communicate with a separate visual stimulus program (`CocoVideoClient.m`)
- Visual stimulus frequencies are mapped to TCP frequencies (3, 13, 23)
- For 'wind' or 'none' conditions, a special value (100) is sent to indicate no visual stimulus

#### Video Recording
- Camera: Dorsal view (labeled as "Coronal (x) Camera")
- Resolution: 640x480 pixels, Y8 format
- Frame rate: 60 fps
- Manual shutter mode with shutter value of 800
- Videos saved as Motion JPEG AVI files
- Camera is hardware-triggered via the DAQ system

### How to read it:

```python
from pymatreader import read_mat
from pathlib import Path

folder_path = Path("/home/heberto/cohen_project/Sample data/Suver Lab/behavioralDataExamples")
file_path =  folder_path / "2024_07_08_E3.mat"
assert file_path.is_file()

mat_data = read_mat(file_path)["data"]
number_of_trials = len(mat_data)
print(f"number_of_trials={number_of_trials}")


non_array_keys = []
array_keys = []

for key, value in mat_data.items():
    first_value = value[0]
    if isinstance(first_value, np.ndarray):
        array_keys.append(key)
        print(f"{key=} {first_value.shape=}")
    elif isinstance(first_value, list):
        print(f"{key=}, len(value)={len(value)}")
        array_keys.append(key)
    else:
        non_array_keys.append(key)
        print(key, first_value)

number_of_trials=15
date None
expnumber 3
condition dark
min_age 3
max_age 5
samplerate 20000
adjust_time 3
record_time 3
fps 60
nframes 2520
trial 1
stimType 1
key='pufferSignal' first_value.shape=(960000,)
key='tachometerSignal' first_value.shape=(960000,)
key='tachometerSignal_smoothed' first_value.shape=(960000,)
```

Note that some of the dates are None, in that case we are adding the date as the current data.
