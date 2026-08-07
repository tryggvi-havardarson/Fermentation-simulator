# Fermentation Simulator

A Python-based batch fermentation simulator for chemical engineering, process modeling, and biochemical reaction simulation.

## Project Goal

The goal of this project is to build a fermentation simulator that models microbial growth and product formation in a batch bioreactor. As the project develops, the simulator will incorporate more realistic biological and engineering phenomena while serving as a tool to strengthen my Python programming and chemical engineering modeling skills.

## Current Status

This project is being developed in multiple versions, with each version focusing on improving both the simulator and my understanding of Python and fermentation modeling.

### Version 1

- Procedural Python implementation.
- Monod growth kinetics.
- Euler numerical integration.
- Biomass, glucose, ethanol, and carbon dioxide prediction.
- Fermentation profile visualization.
- Component mass balance.

### Version 2

- Object-oriented redesign.
- Rosso cardinal temperature model.
- Yeast strain database.
- Improved project structure and code organization.
- Input validation.
- Modular simulation workflow.

## Planned Version 3

- ABV, specific gravity (SG), and Brix calculations.
- Feedstock database with sugar composition.
- Nutrient availability model.
- Target batch calculator (volume, feedstock, and target ABV).
- CO₂ production and fermenter weight calculations.
- Unit conversion tools.
- Process flow diagram.
- Fermentation heat generation model.
- Improved biological realism through growth-phase modeling.


## How to Use

1. Clone or download this repository.
2. Install the required Python packages:
   - NumPy
   - Matplotlib
3. Open `main.py`.
4. Edit the simulation parameters.
5. Run the program.
6. The simulator will display the fermentation profile and print a simulation summary, including a component mass balance.

## Current Assumptions

- Batch fermentation.
- Well-mixed bioreactor.
- Constant reactor temperature.
- Constant reactor volume.

## Current Limitations

- Uses simplified fermentation kinetics.
- Uses proxy yeast parameters.
- Simulation parameters must be edited directly in the source code.
- Does not account for pH, dissolved gases, nutrient limitations, contamination, maintenance metabolism, or changing reactor conditions.