# Fermentation simulator

A Python-based batch fermentation simulator for chemical engineering, process modeling, and biochemical reaction simulation.

## Project Goal

The goal of this project is to build a fermentation simulator that models microbial growth and product formation in a batch bioreactor. As the project develops, the simulator will incorporate more realistic biological and engineering phenomena while serving as a tool to strengthen my Python programming and chemical engineering modelling skills.

## Current Status

**Version:** 1.0.0
This is the first complete release of the fermentation simulator. It models batch fermentation using a simplified kinetic model and provides visualisation of the simulation results. Future versions will expand the model with additional biological and process engineering features.

## How to use

1. Clone or download this repository.
2. Install the required Python packages:
   - NumPy
   - Matplotlib
3. Open `batch_fermentation_simulator.py`.
4. Edit the simulation parameters in the `main()` function if you want to change the fermentation conditions.
5. Run the script.
6. The program will display a graph of the fermentation profile and print a summary of the simulation results.

## Assumptions

- Batch fermentation.
- Well-mixed bioreactor.
- Constant temperature throughout the bioreactor.

## Limitations

- Supports one yeast strain (CBS 8066).
- Uses a simplified kinetic model that does not account for thermal inactivation.
- Simulation parameters must be edited directly in the source code.
- Does not account for pH, dissolved gases, temperature changes during fermentation, contamination, or nutrient limitations.
