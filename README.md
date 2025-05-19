# Circadian Clock-Regulated Vaccine Response Simulation

This repository contains source code and experimental data for simulating circadian clock-regulated vaccine responses.

## Requirements

- Recommended: **Matlab 2021 or later**
- Required Toolboxes:
  - Statistics and Machine Learning Toolbox
  - Navigation Toolbox

## File Specifications

### Source Code

- `demo_fig_3.m`  
  Generates Fig.3 (Simulation results Circadian steady-state migration of immune system, innate immune response, and time-of-day-dependent adaptive immune response)

- `demo_fig4_fig5_fig6.m`  
  Generates Fig.4-6 (Simulation results for three vaccine types)

- `demo_fig8.m`  
  Generates Fig.7 (Dynamic bifurcation diagram of self-accelerating DC homing process)

- `demo_fig7_fig9_figS4.m`  
  Generates Fig.8, Fig.9 and Fig.S4 (Simulated DC homing process under different vaccine types)

- `demo_figS5.m`  
  Sobol sensitivity analysis of circadian-controlled immune parameters on adaptive immunity

- `plot_sobol_indices.m`  
  Script for visualizing Sobol sensitivity indices

- Visualization Utilities:
  - `symlog_2.m` / `symlog_deepseek.m` (Symlog coordinate plotting)
  - `arrowPlot.m` (2D curve annotation with arrows)
  - `arrowPlot3D.m` (3D curve annotation with arrows)

### Parameter Tables

- Core Parameters:
  - `par_clock.csv` (Core circadian clock network parameters from Abo et al. 2021)
  - `circa_immune_list.xlsx` (Circadian-immune regulation parameters)

- Dynamic Process Parameters:
  - `produce_list.xlsx` (Cell activation/differentiation and cytokine production)
  - `decay_list.xlsx` (Component degradation rates)
  - `transition_list.xlsx` (Migration/diffusion coefficients of immune cells and cytokines across compartments)

- Interaction Parameters:
  - `half_maximum_list.xlsx` (Half-maximal effect parameters)
  - `promotion_list.xlsx` (Amplification coefficients)

- Auxiliary Data:
  - `Timing.xlsx` (Temporal profiles of 12 circadian components)

  - `Par_ifam_ada_index.mat` (Weighting coefficients for inflammation/adaptive response indices)

  - `Z.mat` (Fixed-point set for DC-sCCL21 loop)

### Early Model-Experimental Data Fitting

- `055.fig`  
  Simulated LPS-induced cytokine responses (3/6/12 mg/kg) from early models

- `adaptive_simu.fig`  
  Fitting results of adaptive immunity experimental data from early models

    
