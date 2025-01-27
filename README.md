# Polysulfamide Simulations and Analysis

This folder contains code and resources for running and analyzing **polysulfamide** simulations. It includes modifications made to the **atomistic OPLS-AA force field** and tools for the analysis of **self-assembly** in various **homopolymers** and **mixtures**.

## Atomistic Simulations:
- **Atomistic Forcefield Modifications:** Changes made to the OPLS-AA forcefield for polysulfamide simulations.
- **PDB File Generation:** Scripts to generate PDB files for aliphatic backbone polysulfamide.
- **Simulation Setup:** GROMACS commands and settings for running polysulfamide atomistic simulations.
- **Trajectory Analysis:** Tools for analyzing the simulation trajectories and extracting key data.

### How to Use:
1. Clone or download the repository.
2. Generate the PDB file using the Python script in **Create pdb file**, selecting the number of chains and repeat units you want. Then, add the files from **opls_Polysulfamide** to your OPLS folder.
3. Run simulations using the provided **Simulations** folder.
4. Analyze the results using the **Analysis** scripts.
   
## Coarse Grain Simulation Analysis Tools:
Code for analyzing self-assembly behavior, including properties of homopolymers and mixtures.
   
### How to Use:
1. Clone or download the repository.
2. Use the provided analysis scripts to evaluate the CG results.

## Dependencies:
- **Atomistic Simulations:**  
  - GROMACS - [https://www.gromacs.org/](https://www.gromacs.org/)
  
- **Coarse-Grain Simulations:**  
  - LAMMPS - [https://docs.lammps.org/Manual.html](https://docs.lammps.org/Manual.html)
  
- **Simulation Analysis:**  
  - Python with the following libraries:  
    - Numpy  
    - Pandas    
    - MDTraj
    - Seaborn
    - Matplotlib
    
