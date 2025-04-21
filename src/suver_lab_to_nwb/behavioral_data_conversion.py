"""
Suver Lab Behavioral Data Conversion

This script demonstrates how to use the three behavioral data conversion modules:
- windy_steps_conversion.py
- speedy_bars_conversion.py
- coco_conversion.py

Each module converts a specific type of behavioral experiment data to NWB format.
"""

from pathlib import Path
import pandas as pd
from pynwb import read_nwb

from windy_steps_conversion import convert_windy_steps_to_nwb
from speedy_bars_conversion import convert_speedy_bars_to_nwb
from coco_conversion import convert_coco_to_nwb


def main():
    """Run all three behavioral data conversions and display trial tables."""
    # Define common paths
    data_dir = Path("/home/heberto/cohen_project/Sample data/Suver Lab/behavioralDataExamples/")
    output_dir = Path("/home/heberto/cohen_project/Sample data/Suver Lab/nwb")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. WindySteps Conversion
    print("\n" + "="*80)
    print("WINDY STEPS CONVERSION")
    print("="*80)
    print("In this experiment, each fly experienced 10 airflow trials either increasing")
    print("from 0 cm/s to 300 cm/s (5 trials) or decreasing from 300 cm/s to 0 cm/s (5 trials).")
    print("Airflow was presented in 50 cm/s steps over 42 seconds, with 6 seconds of steady")
    print("state airflow for each windspeed.")
    
    windy_steps_file = data_dir / "2024_07_08_E3.mat"
    windy_steps_nwb = convert_windy_steps_to_nwb(
        matlab_data_file_path=windy_steps_file,
        output_dir=output_dir,
        overwrite=True,
        verbose=True,
    )
    
    # Display WindySteps trial table
    nwbfile = read_nwb(windy_steps_nwb)
    trials_df = nwbfile.trials.to_dataframe()
    print("\nWindySteps Trial Table:")
    print(trials_df)
    
    # 2. SpeedyBars Conversion
    print("\n" + "="*80)
    print("SPEEDY BARS CONVERSION")
    print("="*80)
    print("In this experiment, each fly experienced 60 optic flow trials ranging from 0–50 cm/s")
    print("in 5 cm/s increments along with 100 cm/s (5 trials at each optic flow speed,")
    print("pseudo-randomized order). Each trial consisted of an 8 second period of constant optic flow.")
    
    speedy_bars_file = data_dir / "2024_05_29_E3.mat"
    speedy_bars_nwb = convert_speedy_bars_to_nwb(
        matlab_data_file_path=speedy_bars_file,
        output_dir=output_dir,
        overwrite=True,
        verbose=True,
    )
    
    # Display SpeedyBars trial table
    nwbfile = read_nwb(speedy_bars_nwb)
    trials_df = nwbfile.trials.to_dataframe()
    print("\nSpeedyBars Trial Table:")
    print(trials_df)
    
    # 3. Coco Conversion
    print("\n" + "="*80)
    print("COCO CONVERSION")
    print("="*80)
    print("In this experiment, each fly experienced 36 trials at 3 oscillatory frequencies")
    print("(0.3, 1.3, 2.3 Hz, 12 trials at each frequency) with windspeeds ranging from")
    print("100–200 cm/s and optic flow speeds from 5–35 cm/s.")
    
    coco_file = data_dir / "2024_10_28_E4.mat"
    coco_nwb = convert_coco_to_nwb(
        matlab_data_file_path=coco_file,
        output_dir=output_dir,
        overwrite=True,
        verbose=True,
    )
    
    # Display Coco trial table
    nwbfile = read_nwb(coco_nwb)
    trials_df = nwbfile.trials.to_dataframe()
    print("\nCoco Trial Table:")
    print(trials_df)
    
    # Summary
    print("\n" + "="*80)
    print("CONVERSION SUMMARY")
    print("="*80)
    print(f"WindySteps NWB file: {windy_steps_nwb}")
    print(f"SpeedyBars NWB file: {speedy_bars_nwb}")
    print(f"Coco NWB file: {coco_nwb}")
    print("\nAll behavioral data has been successfully converted to NWB format.")


if __name__ == "__main__":
    main()
