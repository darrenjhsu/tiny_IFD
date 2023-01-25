# tinyIFD
Light-weight induced fit docking pipeline through AmberTool's `mdgx.cuda`.

## Requirements

You need the following software.

1. AmberTools
  - `tleap` for parametrizing peptides 
  - `antechamber` for parametrizing the ligand 
  - A modified version of `mdgx.cuda` - for code, see here
  - `cpptraj` for outputting features that the classification models will take
1. Python packages - they should be in a separate environment 
  - `openbabel` for conversion between formats
  - `AutoDock Vina` for docking
  - `rdkit` for various modeling tasks
  - `mdtraj` for easy trajectory manipulation
  - `OpenMM` for energy minimization (a version where we use this as the simulation engine is being developed) 
  - `XGBoost` for classification model
  - `spyrmsd` for symmetry corrected RMSD
  - `oddt` for molecule modeling tasks

## Getting tinyIFD to work for OLCF users

Note that, since we only need cuda when running `mdgx.cuda`, which requires us to be on Summit, the rest of the workflow can be done on Andes.
However, this does require you to install AmberTools on both Andes and Summit.


### Install python packages on Andes

### Install AmberTools on Andes

### Install AmberTools (specifically mdgx) on Summit


## Tutorial with one docking case

### Prepare input files

### Process input files

### Run docking

### Parametrize ligand

### Generate complexes of ligand and truncated protein core

### Energy minimization

### Simulation with mdgx.cuda

### Process trajectories

### The MDR refinement script and data structure

### Predicting docking poses after refinement










