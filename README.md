# Parametrically_Excited_Pre-buckled_Strut_Paper
## A. Overview

This repository contains the data supporting the paper:

📌 Please cite the following paper if you use this repository

Chen, Y., Wadee, M.A., Málaga-Chuquitaype, C. (2026). Nonlinear resonance characterization of parametrically excited pre-buckled elastic struts.
Proceedings of the Royal Society A DOI: 10.1098/rspa.2026.0072

## B. Files in the repository

The dataset follows a consistent naming convention to encode both the models and their corresponding results.

1. Figures: This folder contains PDF versions of all figures presented in the main article and the supplementary material, using identical notation to that in the
manuscript.

2. Analytical data (MATLAB, MATCONT):
   
   (a) Function and script: This folder contains .m files defining the dynamical systems within the MATCONT framework, referred to as odefiles. Each odefile includes the governing equations of the dynamical system, the associated ini tialisation parameters, and the symbolic derivatives required by MATCONT. The odefile alone is sufficient for use with the command-line version of MATCONT. For the graphical user interface (GUI) version, both the corresponding odefile and a .mat file with the same name must be placed in the system folder of MATCONT. The naming convention for the odefiles is as follows:
   
   i. BmDyn1m3t: first-mode discretised, cubic-order parametrically excited system (the primary model used in the article);
   ii. BmDyn13m3t: first- and third-mode discretised, cubic-order parametrically excited system;
   iii. BmDyn1m3t: first-mode discretised, seventh-order parametrically excited system;
   iv. BmDynD1m3t: first-mode discretised, cubic-order directly excited system.
   
  (b) Results: This folder contains the results corresponding to the models listed above.

3. Finite element data (ABAQUS):
   
(a) Script: This folder contains the Python scripts used to generate the ABAQUS input files for an axially excited strut, together with an example input file.

(b) Results: This folder contains multiple sets of steady-state results obtained from the finite element model, stored in .mat format.

4. README.pdf: file describing the nomenclature used for the files in this repository.

5. Suplementary_Material.pdf:
   S1 Directly excited case: Nonlinear resonance characteristics.
   S2 Characteristics of the PPR peak: FE verification.
   S3 Poincar´e map and bifurcation diagram.
   S4 Significance of axial inertia nonlinearity. 

## C. Repository DOI

You can cite this repository as follows:  Chen, Y., Wadee, M.A., Málaga-Chuquitaype, C. (2026). Shared Data and Supplementary Material of 'Nonlinear Resonance Characterization of Parametrically Excited Pre-buckled Elastic Struts'. Zenodo, DOI: 10.5281/zenodo.19701820
