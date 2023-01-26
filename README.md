# tinyIFD
Light-weight induced fit docking pipeline through AmberTool's `mdgx.cuda`.

## Requirements

You need the following software.

1. AmberTools
    - `tleap` for parametrizing peptides 
    - `antechamber` for parametrizing the ligand 
    - A modified version of `mdgx.cuda` - for code, see [here](https://github.com/darrenjhsu/mdgx_mod/tree/fix_atom) and below
    - `cpptraj` for outputting features that the classification models will take
1. Python environment manager (conda highly recommended)
1. Python packages - they should be in a separate environment
    - `openbabel` for conversion between formats
    - `AutoDock Vina` for docking
    - `rdkit` for various modeling tasks
    - `mdtraj` for easy trajectory manipulation
    - `OpenMM` for energy minimization (a version where we use this as the simulation engine is being developed) 
    - `XGBoost` for classification model
    - `spyrmsd` for symmetry corrected RMSD
    - `oddt` for molecule modeling tasks
    - `parmed` for molecule modeling tasks
    - `numba` for parallel processing
  

## Getting tinyIFD to work for OLCF users

Note that, since we only need cuda when running `mdgx.cuda`, which requires us to be on Summit, 
the rest of the workflow can be done on Andes.
However, this does require you to install AmberTools on both Andes and Summit.
For the python packages, all of them can be on Andes.
A python interpreter on Summit is necessary, but no fancy packages are required.

### Install python packages on Andes

First, log on to Andes.

```bash
ssh user@andes.olcf.ornl.gov
```

Load modules that will be used by AmberTools

```bash
module load gcc/9.3.0 cmake boost \
            netcdf-c netcdf-cxx netcdf-fortran parallel-netcdf \
            openblas netlib-lapack fftw
```

Load python and conda

```bash
module load python
conda init bash
source ~/.bashrc
```

Now create an environment for tinyIFD.

```bash
conda create -p $PWD/myenv python=3.8
conda activate $PWD/myenv
```

Install python packages

```bash
conda install numpy=1.23 openbabel vina rdkit mdtraj openmm xgboost spyrmsd oddt parmed=3.4.3 numba
```

Note that, `numpy` 1.23 is required for `np.int` in `vina`, and `parmed` 3.4.3 is required to *not* write `CONECT` records when outputting pdb files.

### Install AmberTools on Andes

Now, let's install AmberTools on Andes. 

1. Get AmberTools source code (in my case it is called `AmberTools20.tar.bz2` and upon decompression it gives `amber20_src/`) from the website.
1. Go to `amber20_src/AmberTools/src` and clone the modified mdgx code:

    ```bash
    cd amber20_src/AmberTools/src
    git clone git@github.com:darrenjhsu/mdgx_mod.git
    ```

1. Backup the original mdgx code and swap in this modified version

    ```bash
    mv mdgx mdgx_bkp
    mv mdgx_mod mdgx
    cd mdgx
    git checkout fix_atom
    ```

1. Use `ccmake` to install AmberTools from `amber20_src`

    ```bash
    # In amber20_src
    mkdir build_andes
    cd build_andes
    ccmake ..
    ```

In the interactive panel, first do [c] to configure, then change the following:

```
CMAKE_INSTALL_PREFIX            */path/to/amber20_src/build_andes
COMPILER                        *GNU
```

Do [c] again, then 

```
DOWNLOAD_MINICONDA      OFF
OPENMP                  ON
```

Do [c] again, then

```
BUILD_PYTHON            OFF
```

Do [c] again, and then [g] to generate configs. After this step, do

```
make -j 16 | tee build.log
make install | tee -a build.log
```

This should install the AmberTools on Andes. Remember to include in your `PATH` the `bin/` folder.

```
export PATH=$PWD/bin:$PATH
```

### Install AmberTools (specifically mdgx) on Summit

Please follow [this tutorial](https://github.com/darrenjhsu/mdgx_on_Summit) to do so.

Also, it is helpful to have at least a python interpreter on Summit:

```bash
module load python
```

## Tutorial with one docking case

In this tutorial we are simply going to dock and refine the ligand 
from the SARS-CoV-2 main protease 7QBB to the structure of 5R84, 
a typical "cross-docking" job. Because the active site in 5R84 accommodates to 
its own ligand, simple docking (such as by Vina) may not yield a correct structure.

Through this tutorial, we are going to learn how to move data through the tinyIFD pipeline.

### File structure

As you can see, before we run anything, there are `01_Workflow/` and `02_Input/` folders. 
The folder `01_Workflow/` holds scripts that read input files from `02_Input/` folder,
which we prepare, and generate scripts to process them.

Inside `02_Input/`, we have three folders `Receptors/`, `Ligands/`, and `Ref_Receptors/`.
Within each, there is a single folder `7QBB-V1B-to-5R84/`, meaning "we dock the ligand V1B from 7QBB to 5R84".

`Receptors/7QBB-V1B-to-5R84/` contains the "template apo receptor" which is the one **to be docked to**, or 5R84.
`Ligands/7QBB-V1B-to-5R84/` contains the "target ligand" which is the ligand from 7QBB.
It also contains the template ligand, which is the ligand from 5R84.
This is used in some of the features downstream in the pipeline.
`Ref_Receptors/7QBB-V1B-to-5R84/` contains the "target holo receptor" which is the one **where the ligand is from**, or 7QBB.

Finally, in `02_Input/`, a `job_description.csv` file denotes all the tasks to be done, file locations, as well as the docking centers.
The filename is hardcoded for now. 


### Prepare input files

I have prepared this one, but in general you want to align the target PDB to template PDB, and use the center of your template receptor's ligand as your docking center. A separate tutorial is in development.

### Process input files

Now let's get some work done. Head into `01_Workflow/` and do the following

```bash
python 00_Parse_job_description.py -h   # Prints parameters that you can supply
python 00_Parse_job_description.py --config 00_config  # Reads 00_config for config parameters; then reads 02_Input/job_description.csv and check that all file exists
```

The `00_config` file just denotes resources required and some stretagies.

```
dockingMode=rigid           # Do rigid docking; flexible docking is under development
parallelGrid=32             # Calculate grids. This version of tinyIFD doesn't actually do it. Just keep it as is
parallelDock=1              # How many Andes nodes you need to dock things (e.g. 25 nodes for 100 cases)
assembleMode=protein_atoms  # Assemble mode for the tiny system - protein_atoms, atoms, or residues are possible choices
assembleTarget=800          # Try to contain about 800 protein atoms, then add the ligand
parallelAntechamber=1       # Andes nodes you need for parametrizing ligand (similar to parallelDock)
parallelPrepareComplex=32   # Keep this at 32
parallelEM=1                # Andes nodes for OpenMM energy minimization
parallelcpptraj=1           # Andes nodes
parallelMDR=1               # Andes nodes for reducing trajectories
parallelXGBoost=1           # Andes nodes for creating XGBoost models
```

Those that are 1 are parallelized at the node level (each task uses 32 cores), 
while those that are 32 are parallelized at the core level (each task uses 1 core).

Now go back to the root, you will see `03_Gridmaps/`, `04_Docking/`, `05_Refinement/`, and `06_Analysis/` folders.
We will go into each folder and run something.

### Prepare protein (PDB -> PDBQT)

In the past versions of tinyIFD we used to use AutoDock 4, and had to generate grids.
We have since switched to the python AutoDock Vina and it makes things much easier.
However, we still need the `.pdbqt` protein coordinates.

```bash
cd 03_Gridmaps/script
sh prepDock0.sh        # Cleans .pdb files and generate .pdbqt files
```

### Run docking

Then, we perform actual docking in `04_Docking/`

```bash
cd ../../
cd 04_Docking/script
sh dock0.sh            # Docks and prepares the protein core
sh check_docking.sh    # Check docking result, RMSD in Angstroms. This is only meaningful when the provided ligand is the crystal structure
```

### Parametrize ligand

After docking and before running MD, we need to parametrize ligand through antechamber.

```bash
cd ../../
cd 05_Refinement/script
sh antechamber0.sh          # Or in the case of multiple tasks, do sbatch andes_antechamber.sh
```

### Generate complexes of ligand and truncated protein core

Then we prepare complexes. Also, to use OpenMM's syntax for freezing atoms, we run a separate script.

```bash
sh prepareComplex0.sh       # Or sbatch andes_prepareComplex.sh, and skip next line
sh prepareComplex_openmm0.sh
sh check_atom_num.sh        # Check number of atoms in each system
```

### Energy minimization

Use OpenMM to minimize energy of the system

```bash
sh energyMinimization0.sh   # Or sbatch andes_EM.sh
```

### Simulation with mdgx.cuda

Now the fun begins! The `planEM.sh` is essentially a planner looking at how many systems you have.
It will ask you how many GPUs you want to use, but for this phase we just do 1.

```bash
sh planEM.sh                # Enter 1 for using 1 GPU
sh EM_script.sh             # If you look inside, it just creates folders, but this script is important! Otherwise mdgx.cuda will segfault
```

And then, the majority of compute actually happens here. 
I have lowered the `N-rep` in the `MD` stage in the `simulation_config.json` to 4
and simulation to 2 ns (1 M steps)
so that all compute fits in 1 GPU and finishes quickly. 
Ideally you want 20 replicas per pose and 8 ns simulations. 

```bash
vim ../../01_Workflow/utilites/simulation_config.json # Read, and see what kind of parameters you can adjust
sh planMD_openmm.sh         # Enter 1 for using 1 GPU
```

This generates `mdgxGPU_MD_0.in`, which is the input to `mdgx.cuda`, as well as `Ssubmit_MD_0.sh`, 
which is the job submission file for summit. 
Let's go to summit and actually do it.

```bash
# On summit
sh MD_local.sh
# bsub -q debug Ssubmit_MD_0.sh # If you need more than 1 GPU, then do this
```

This will take about 25 minutes to complete on 1 GPU


### Process trajectories

With the simulation done, we can check if there are divergent simulations that failed.
It is alright to have some (like 4 - 5 for a total of 400).

```bash
$ python check_sims.py

There are 1 IFD cases.

Ligand: 7QBB-V1B-to-5R84 | 80 simulations | 80 finished | 80 succeeded | 0 failed (0.0 %)
Total simulations: 80
Total failed simulations: 0
Total success simulations: 80
Success percentage: 1.0
```

Now, let's process the trajectories on Andes.
First, generate the `cpptraj` scripts for each of the simulations

```bash
sh gen_cpptraj_script.sh 
```

Then, depending on how many simulations you have, either do
```bash
# For not a lot of simulations
sh run_cpptraj_script.sh
```

or do 

```bash
# For A LOT OF simulations
sbatch run_cpptraj_script_andes.sh
```

### The MDR refinement script and data structure

With the data written by `cpptraj`, we then digest the entire set of simulations with the MDR data structure:

```bash
cd ../../06_Analysis/script
sh MDR_analysis_all.sh  # Sequentially go through each docking job
# or 
# sbatch MDR_analysis_andes.sh  # Go through docking jobs in batches
```

This would take some time, usually 5 to 6 mins per round.

### Predicting docking poses after refinement

Under construction.







