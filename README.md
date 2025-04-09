## Motor Cortex Activity Analysis Using Raster and PETH Plots
### Description
This repository contains MATLAB code and analysis to visualize and interpret neuronal activity from a motor cortex array in response to various finger movements. Using Raster plots and Peri-Event Time Histograms (PETHs), this project explores how specific neurons encode distinct finger movements.
### Data Overview
The analysis uses electrophysiological recordings from a 10x10 microelectrode array implanted in a macaque monkey’s motor cortex. The neural activity was recorded during specific finger movements: thumb, index, middle, and combinations thereof.
### Analysis Techniques
- Raster Plots: Visualization of neuron firing aligned to behavioral events.
- PETH Plots: Calculation and plotting of neuron firing rates around specific events (cue or actual movement).
- Decoding Criteria: Established criteria based on firing rates to decode finger movements from neuronal activity.
### Repository Structure
- data/ - Contains neural and behavioral data (.mat files).
- scripts/ - MATLAB scripts performing Raster and PETH analysis.
- figures/ - Generated visualizations showing neuronal activity aligned to finger movements.
### Dependencies
- MATLAB environment
- Basic MATLAB toolboxes

### Contributions
Contributions and improvements are welcome. Please fork the repository, submit pull requests, or raise issues for discussion.
