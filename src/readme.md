# Monte Carlo Simulation of Photon Transport in Planetary Atmospheres

## Overview
This project implements a Monte Carlo Radiative Transfer (MCRT) simulation to study the scattering of photons in planetary atmospheres, particularly relevant for bodies like Titan. The method traces individual photons through interaction processes, including absorption, scattering, and surface reflection, using probabilistic models.

## Files
- `monto curlo .docx`: Main report with methodology and results.
- `Monte-Carlo-Project-Overview.pptx`: Presentation slides summarizing the approach.
- `three mcm.pdf`: Supplemental visuals including flowcharts and result interpretation.
- `output_path.jpeg`: Visual results showing photon paths and intensity distributions.

## Requirements
- Python 3.x (with NumPy and Matplotlib)
- Alternatively, Fortran/C++ can be used with proper random number generation and ray-tracing logic.

## How to Run
1. Initialize the atmospheric model (layer properties, albedo, phase function).
2. Launch photons from the top of the atmosphere.
3. Track each photon's path using probabilistic rules until absorption, escape, or surface interaction.
4. Record statistics like intensity distribution, angular dependence, etc.
5. After succesfully generate dat file, it will have 16 file's.
6. import date file's in notebook and run the jupyter notebook.
7. Plot results using Matplotlib.

## Output
- 2D/3D plots of photon paths(plote as you need)
- Histograms of angular distribution
- Comparison against deterministic models like DISORT

## Note:
- We have ploted result for UV wave you may try with IR or other as you needed
- you may also run exe file in case of you unable to run program
