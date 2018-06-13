# Master simulator

This is a framework for simulating phylogenetic trees in MASTER using BEAST2 and then inferring parameters on the simulated trees using different models. 

The framework is implemented in Python and relies on template XML files for the MASTER simulations and BEAST inferences, as well as config files in YAML for setting model and inference parameters.

The scripts are compatible with BEAST 2.5.0 and MASTER 6.0.0


## File formats

- **Config files:** Standard YAML files, with variable names and values. Comments are allowed. Some special keywords are recognised that do not appear in the template XML files:
	- `template`: The template XML file to use 
	- `name`: The name of the output XML file (will be `name.xml`)
	- `outputpath`: The directory to store the output XML file in

- **Template files:** Standard BEAST2 XML files, with hooks to be replaced specified as e.g. `{$var}`, which will be replaced by the value associated with `var` in the config file.


## Scripts

### 1.) MakeMasterXML.py
Creates a MASTER XML file. Input parameters as well as the template XML file should be specified in a YAML config file. If an input directory is supplied, the script will be run on every file in the directory that matches the pattern for a config file (default = `*.cfg`). 

Some parameter conversions for birth-death-sampling type models (as implemented in BEAST2) are hardcoded into the script, so only one parameterisation needs to be included in the config file. In particular to add in both parameterisations for the birth-death skyline model (as infection/recovery/sampling rates and as reproductiveNumber/becomeUninfectious/samplingProportion). 


### 2.) MakeBeastXML.py
Creates a BEAST2 XML file. Input parameters for the template XML file should be specified in a YAML config file. To make it possible to use the same config file as the simulations the template file, output directory and filename may be specified as input options. If an input directory is supplied, the script will be run on every file in the directory that matches the pattern for a config file (default = `*.cfg`). 

The script will look for the output trees of the simulation and iterate through all trees in the simulation output, creating a BEAST2 XML file for each tree. The sampling times of tips in the tree are read from the output file and the tMRCA is set to time = 0.

The script also produces bash scripts for running inferences for each replicate of a simulation, as well as a script for the ETH Zürich cluster (Euler). 


### 3.) CheckBEAST.py
Need to check if this works the way it is supposed to (summarise results from replicate simulations and automatically check convergence).


### 4.) plotSimulations.R
Plots the result of the MASTER simulations in R. This file is currently hard-coded for the 3 types of simulations in the project, but can be easily extended.


## Example workflow

```bash

# Create MASTER XML files
python MakeMasterXML.py -i ../config/

# Run MASTER simulations
cd ../results/simulations/
for i in `ls *.xml`; do ~/Documents/Packages/BEASTv2.5.0/bin/beast $i; done
cd -

# Plot simulation results in R
# (run plotSimulations.R)

# Create BEAST XML files (birth-death skyline)
python MakeBeastXML.py -i ../config/ -o ../results/inferences/ -x ../templates/inferences/template.bdsky.fixtree.xml -n bdsky.fixtree

# Create BEAST XML files (Bayesian skyline plot)
python MakeBeastXML.py -i ../config/ -o ../results/inferences/ -x ../templates/inferences/template.bsp.fixtree.xml -n bsp.fixtree

# Run inferences


# Check convergence and create summary statistics across replicates

```


## Disclaimer

These are scripts from an unfinished project that I thought should probably be online and would be of use to other people. They are provided as is and are not guaranteed to be correct.