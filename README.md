# SCP-dielectric-model
MATLAB implementation of a physically based soil complex permittivity (SCP) model, including PBm solver, DCC parameterization, and corresponding LUTs.


This work develops a soil dielectric modeling framework that consists of three main components, arranged in a sequential workflow:
(i) a Dynamic Cole–Cole (DCC) model to describe the dielectric response of saline solutions,
(ii) a Poisson–Bikerman (PBm) model to explicitly account for bound-water and electric double layer effects, and
(iii) an SCP model to characterize the overall complex permittivity of saline soils.

Accordingly, the core code for the saline-solution DCC dielectric model, the core implementations and data for solving the PBm model via a boundary value problem (BVP) approach, and the core code and data for generating SCP-based dielectric properties are organized into three subdirectories: DCC, PBm_bvp, and SCP_LUTs, respectively.


## Basic Workflow

1. Fit the DCC model to saline-solution measurements using scripts in `DCC/`.
2. Solve the PBm model to obtain ion distributions and local dielectric profiles (`PBm_bvp/`).
3. Generate SCP-based dielectric properties and LUTs for soil samples (`SCP_LUTs/`).

- `DCC/`  Implementation of the DCC model for saline solutions.  

  The main fitting routine is provided in `fit_cc_quadratic_c.m`. As the original experimental data are not publicly available, the fitted outputs are included in `best_quad_cc_equations.txt`, namely the optimized ε_∞, and the analytical expressions of Δε, τ, and β as the functions of concentration c.
  

- `PBm_bvp/`  Numerical solvers for the PBm model using bvp4c methods.

  This folder contains two main scripts: `PBm_main.m`, which serves as the main driver, and `new_solveEDL_sigma_iter_iterative.m`, which implements the core numerical solver for the PBm model.

  Using the DY soil sample as an example, running the main script produces the bound-water layer thickness d corresponding to different initial concentrations c₀, which are saved in `all_dthick_nm.txt`. The `results/` directory further stores the ion distribution profiles within the bound-water layer, as well as the associated electric potential and electric field intensity for different combinations of soil salinity (SS) and soil moisture (SM) (see in-code comments for details).

  To adapt the model to other soil types, the initial soil-solution concentration c₀ should first be computed for different SS and SM conditions and saved as `all_mol.txt`, and the clay content parameter in the main script should then be updated to the target soil value before running the model.

  

- `SCP_LUTs/`  Scripts and lookup tables (LUTs) for generating SCP-based dielectric properties for various soil samples.
  
  Running `generate_SCP_LUT_all_samples.m` directly will generate SCP-based lookup tables using the embedded parameter settings for the six soil samples considered in this study. The script produces the complex permittivity (ε) LUTs for both the L band (1.4 GHz) and the C band (5.4 GHz) for all six samples, which are saved in the `LUT/` subdirectory.

  If permittivity LUTs at other frequencies within the 1–6 GHz range are required, the corresponding frequency values can be modified directly in the script.
