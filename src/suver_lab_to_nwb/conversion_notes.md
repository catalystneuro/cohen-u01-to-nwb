
## Data organization

```
[4.0K]  .
├── [ 74M]  2024_10_24_E2.mat
├── [ 28K]  DeepLabCut
├── [4.0K]  seal tests
└── [ 20K]  videos
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

The content of one of the fiels looks like this

```python
date None
expNumber 6
trial 1
notes fly01
SAMPLERATE 10000
TRIAL_TIME_LIGHT 4
PRE_TRIAL_TIME 2
POST_TRIAL_TIME 4
numSecOut 10
genotype 24C06
fps 60
nframes 600
key='stimTiming' first_value.shape=(100000,)
key='micLeft' first_value.shape=(100000,)
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
* For the patch clamp, why does the data says 1,000 frames but the length of the data is `100,000 
* samplerate: 10000, fps: 100
* is dep days post eclosing (e.g. 6-8dpe) or days post patching?
* What is scale current?
* What is scale voltage?
* What does the tachometer measures or represents? It is usually used to measure the speed of a rotating object, such as a motor or other machine.
* What about the puffer?
* Description and units of `Vm`, `Im`, `tachometer`, `puffer`, `filteredVm`

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

