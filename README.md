# cohen-u01-to-nwb

NWB conversion scripts for multiple labs' data to the [Neurodata Without Borders](https://nwb-overview.readthedocs.io/) data format. This repository contains conversion tools for data from the following laboratories:

- Marie Suver Lab (University of Washington)
- Itai Cohen Lab (Cornell University)
- Sung Soo Kim Lab (University of California, Santa Barbara)
- Brad Dickerson Lab (University of North Carolina at Chapel Hill)
- Jessica Fox Lab (Case Western Reserve University)

## Installation from GitHub

To install the package directly from GitHub, you will need to use `git` ([installation instructions](https://github.com/git-guides/install-git)). We also recommend the installation of `conda` ([installation instructions](https://docs.conda.io/en/latest/miniconda.html)) as it contains all the required machinery in a single and simple install.

From a terminal (note that conda should install one in your system) you can do the following:

```bash
git clone https://github.com/catalystneuro/cohen-u01-to-nwb
cd cohen-u01-to-nwb
conda env create --file make_env.yml
conda activate cohen-u01-to-nwb-env
```

This creates a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) which isolates the conversion code from your system libraries. We recommend that you run all your conversion related tasks and analysis from the created environment in order to minimize issues related to package dependencies.

Alternatively, if you want to avoid conda altogether (for example if you use another virtual environment tool) you can install the repository with the following commands using only pip:

```bash
git clone https://github.com/catalystneuro/cohen-u01-to-nwb
cd cohen-u01-to-nwb
pip install -e .
```

Note: both of the methods above install the repository in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs).

## Running Conversions

Each lab has its own conversion script(s) that can be run to convert data to NWB format. These scripts serve as blueprints for how common experimental setups in each lab can be converted to the standardized NWB format. They are designed to be readily modified and expanded to meet evolving laboratory requirements. The scripts not only convert primary data streams but also systematically capture experimental, subject, and device metadata, providing practical examples that can be adapted for similar future workflows.

> **Note:** Concrete examples of how to call each conversion function with real data can be found at the bottom of the respective conversion scripts under the `if __name__ == "__main__":` directive. These examples demonstrate typical usage patterns with actual file paths and parameters.

Below are links to the conversion scripts and instructions on how to run them.

### Cohen Lab Conversions (Itai Cohen Lab, Cornell University)

The Cohen Lab has two main types of data conversions:

#### 1. Free Flight Optogenetics Conversion

This conversion processes free flight optogenetics experiments with high-speed video recordings of Drosophila during optogenetic stimulation. It captures:

- Wing and body kinematics (angles, positions)
- Optogenetic stimulation parameters
- High-speed video data (8000 fps)
- Subject metadata (genotype, age, sex)

Conversion script: [free_flight_optogenetics_conversion.py](src/cohen_lab_to_nwb/free_flight_optogenetics_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from cohen_lab_to_nwb.free_flight_optogenetics_conversion import convert_stimuli_experiment_to_nwb
```

#### 2. Confocal Microscopy Conversion

This conversion processes Zeiss confocal microscopy data (.czi files) from immunohistochemistry experiments. It captures:

- Z-stack image data
- Multiple fluorescence channels
- Microscope metadata
- Image dimensions and resolution

Conversion script: [confocal_conversion.py](src/cohen_lab_to_nwb/confocal_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from cohen_lab_to_nwb.confocal_conversion import convert_confocal_to_nwb
```

### Dickerson Lab Conversion (Brad Dickerson Lab, UNC Chapel Hill)

The Dickerson Lab conversion processes Thor imaging data and behavior data from experiments studying neural activity in Drosophila. It captures:

- Two-photon calcium imaging data (multiple channels)
- Behavioral measurements
- Stimulus presentation information
- Experimental metadata

Conversion script: [payel_conversion.py](src/dickerson_lab_to_nwb/payel_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from dickerson_lab_to_nwb.payel_conversion import convert_session
```

### Fox Lab Conversion (Jessica Fox Lab, Case Western Reserve University)

The Fox Lab conversion processes trial data from experiments studying insect flight mechanics and sensory integration. It captures:

- Multiple synchronized camera recordings (side, top, and haltere views)
- DeepLabCut pose estimation data for tracking body parts
- Behavioral measurements and stimulus presentations
- Trial and experimental condition metadata

Conversion script: [streets_conversion.py](src/fox_lab_to_nwb/streets_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from fox_lab_to_nwb.streets_conversion import run_trial_conversion
```

### Kim Lab Conversion (Sung Soo Kim Lab, UC Santa Barbara)

The Kim Lab conversion processes experimental sessions studying sensorimotor integration in Drosophila. It captures:

- Two-photon calcium imaging data
- Behavioral measurements (wingbeat amplitude, flight dynamics)
- Visual stimuli presentation data
- ROI fluorescence traces and coordinates
- Video recordings of fly behavior

Conversion script: [kim_conversion.py](src/kim_lab_to_nwb/kim_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from kim_lab_to_nwb.kim_conversion import convert_session_to_nwb
```

### Suver Lab Conversion (Marie Suver Lab, University of Washington) - Work in Progress

The Suver Lab conversion processes session data with optogenetic stimulation and video recordings. This conversion is currently a work in progress.

Conversion script: [suver_conversion.py](src/suver_lab_to_nwb/suver_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from suver_lab_to_nwb.suver_conversion import convert_session
```

## Repository Structure

The repository is organized by lab, with each lab having its own directory in the `src` directory:

```
cohen-u01-to-nwb/
├── LICENSE                         # Project license file
├── make_env.yml                    # Conda environment specification
├── pyproject.toml                  # Python package configuration
├── README.md                       # This documentation file
└── src/
    ├── cohen_lab_to_nwb/           # Itai Cohen Lab (Cornell University)
    │   ├── confocal_conversion.py
    │   ├── conversion_notes.md
    │   ├── custom_metadata.yml
    │   ├── free_flight_optogenetics_conversion.py
    │   ├── phantom_video_interface.py
    │   └── zeiss_confocal_interface.py
    ├── dickerson_lab_to_nwb/       # Brad Dickerson Lab (UNC Chapel Hill)
    │   ├── behavior_interface.py
    │   ├── conversion_notes.md
    │   ├── payel_conversion.py
    │   └── thor_interface.py
    ├── fox_lab_to_nwb/             # Jessica Fox Lab (Case Western Reserve)
    │   ├── behavior.py
    │   ├── camera_utilites.py
    │   ├── conversion_notes.md
    │   ├── metadata.yaml
    │   └── streets_conversion.py
    ├── kim_lab_to_nwb/             # Sung Soo Kim Lab (UC Santa Barbara)
    │   ├── behavior.py
    │   ├── conversion_notes.md
    │   ├── kim_conversion.py
    │   ├── metadata.yaml
    │   ├── ophys.py
    │   ├── stimuli.py
    │   ├── trials.py
    │   └── utils.py
    └── suver_lab_to_nwb/           # Marie Suver Lab (University of Washington)
    │   ├── conversion_notes.md
    │   └── suver_conversion.py
```

Each lab directory contains:
- Conversion scripts specific to that lab's data.
- Custom interfaces for handling transformation of specific data types / formats to NWB.
- Metadata files (YAML, JSON) for configuring the conversion.
