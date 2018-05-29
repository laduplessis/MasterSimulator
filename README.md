# Master simulator

This is a framework for simulating phylogenetic trees in MASTER using BEAST2 and then inferring parameters on the simulated trees using various models. The framework is implemented in Python and relies on template XML files for the MASTER simulations and BEAST inferences, as well as config files in YAML for setting model and inference parameters.

## Pipeline steps (fixed tree inferences)

1. `MakeMasterXML.py`: This uses the template MASTER XML file along with a YAML config file to create an XML file that will simulate trees using MASTER when run in BEAST2.
2. `MakeBEASTXML.py`: This uses a template BEAST XML file along with a YAML config file and a trees file (produced through MASTER, or any other tool) to create a set of XMLs (one per tree) that can be run in BEAST to infer parameters.
