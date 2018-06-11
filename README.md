# Master simulator

This is a framework for simulating phylogenetic trees in MASTER using BEAST2 and then inferring parameters on the simulated trees using variants of the birth-death skyline model. 

The framework is implemented in Python and relies on template XML files for the MASTER simulations and BEAST inferences, as well as config files in YAML for setting model and inference parameters.


## File formats

- **Config files:** Standard YAML files, with variable names and values. Comments are allowed. Some special keywords are recognised that do not appear in the template XML files:
	- `template`: The template XML file to use 
	- `name`: The name of the output XML file (will be `name.xml`)
	- `outputpath`: The directory to store the output XML file in

- **Template files:** Standard BEAST2 XML files, with hooks to be replaced specified as e.g. `{$var}`, which will be replaced by the value associated with `var` in the config file.


## Scripts

### 1.) MakeMasterXML.py
Creates a MASTER XML file. Input parameters as well as the template XML file should be specified in a YAML config file. If an input directory is supplied, the script will be run on every file in the directory that matches the pattern for a config file (default = `*.cfg`). 

Some parameter conversions are hardcoded into the script, so only one parameterisation needs to be included in the config file. In particular to add in both parameterisations for the birth-death skyline model (as infection/recovery/sampling rates and as reproductiveNumber/becomeUninfectious/samplingProportion). 


### 2.) MakeBEASTXML.py
Creates a BEAST2 XML file


 a template BEAST XML file along with a YAML config file and a trees file (produced through MASTER, or any other tool) to create a set of XMLs (one per tree) that can be run in BEAST to infer parameters.
