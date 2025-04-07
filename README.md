# cohen-u01-to-nwb

NWB conversion scripts for multiple labs' data to the [Neurodata Without Borders](https://nwb-overview.readthedocs.io/) data format. This repository contains conversion tools for data from the following laboratories:

- Marie Suver Lab (Vanderbilt University)
- Itai Cohen Lab (Cornell University)
- Sung Soo Kim Lab (University of California, Santa Barbara)
- Brad Dickerson Lab (Princenton University)
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
- Behavioral video recordings
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

- Two-photon calcium imaging data with GCaMP and tdTomato indicators acquired using the ThorImageLS system
- Detailed wing kinematics:
  - Left and right wing beat amplitude
  - Left minus right wing beat amplitude (differential signal)
- Visual stimulus tracking in X and Y directions
- Experimental metadata and imaging parameters

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
- Detailed wing kinematics:
  - Wing beat amplitude for left and right wings
  - Wing beat frequency
  - Wing voltage recordings showing wing position throughout the stroke
- Synchronization of video, pose estimation and behavioral data.

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
- Synchronization of modalities (imaging, behavior, and stimuli)
- Experimental metadata (genotype, age, sex)

Conversion script: [kim_conversion.py](src/kim_lab_to_nwb/kim_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from kim_lab_to_nwb.kim_conversion import convert_session_to_nwb
```

### Suver Lab Conversion (Marie Suver Lab, Vanderbilt University)

The Suver Lab conversion processes session data from experiments studying Drosophila flight control and sensory integration. It captures:

- In-vivo whole-cell patch clamp recordings (current and voltage)
- Filtered membrane potential recordings
- Tachometer data for detecting wing flapping/flight
- Puffer stimulus data (air puff sensory stimulus)
- Video recordings from multiple angles (lateral fly left, lateral fly right, lateral ventral)
- Pose estimation with DeepLabCut
- Seal test data with calculated metrics (input resistance, access resistance)
- Subject metadata (genotype, age, sex)

Conversion script: [suver_conversion.py](src/suver_lab_to_nwb/suver_conversion.py)

To run this conversion, you can import the conversion function and call it with the appropriate parameters:

```python
# Import the conversion function
from suver_lab_to_nwb.suver_conversion import convert_session_to_nwb
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
    │   ├── confocal_conversion.py  # Converts Zeiss confocal microscopy data
    │   ├── conversion_notes.md     # Documentation of data structure
    │   ├── custom_metadata.yml     # Lab-specific metadata configuration
    │   ├── free_flight_optogenetics_conversion.py  # Converts free flight experiments
    │   └── zeiss_confocal_interface.py  # Interface for Zeiss confocal microscopes
    ├── dickerson_lab_to_nwb/       # Brad Dickerson Lab (UNC Chapel Hill)
    │   ├── behavior_interface.py   # Interface for behavioral data
    │   ├── conversion_notes.md     # Documentation of data structure
    │   ├── payel_conversion.py     # Main conversion script for Dickerson lab data
    │   └── thor_interface.py       # Interface for Thor imaging systems
    ├── fox_lab_to_nwb/             # Jessica Fox Lab (Case Western Reserve)
    │   ├── behavior.py             # Behavioral data processing
    │   ├── camera_utilites.py      # Utilities for camera data processing
    │   ├── conversion_notes.md     # Documentation of data structure
    │   ├── metadata.yaml           # Lab-specific metadata configuration
    │   └── streets_conversion.py   # Main conversion script for Fox lab data
    ├── kim_lab_to_nwb/             # Sung Soo Kim Lab (UC Santa Barbara)
    │   ├── behavior.py             # Behavioral data processing
    │   ├── conversion_notes.md     # Documentation of data structure
    │   ├── kim_conversion.py       # Main conversion script for Kim lab data
    │   ├── metadata.yaml           # Lab-specific metadata configuration
    │   ├── ophys.py                # Optical physiology data processing
    │   ├── stimuli.py              # Visual stimuli data processing
    │   ├── trials.py               # Trial data processing
    │   └── utils.py                # Utility functions for Kim lab conversions
    └── suver_lab_to_nwb/           # Marie Suver Lab (University of Washington)
    │   ├── conversion_notes.md     # Documentation of data structure
    │   └── suver_conversion.py     # Main conversion script for Suver lab data
```

Each lab directory contains:
- Conversion scripts specific to that lab's data.
- Custom interfaces for handling transformation of specific data types / formats to NWB.
- Metadata files (YAML, JSON) for configuring the conversion.
