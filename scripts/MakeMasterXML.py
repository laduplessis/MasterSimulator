import os, sys, yaml
from fnmatch import fnmatch
from optparse import OptionParser

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input files [required]")

parser.add_option("-c","--config",
                  dest = "config",
                  default = "*.cfg",
                  metavar = "",
                  help = "Pattern to match for config files [default = %default]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"

################################################################################################################################  

def makeXMLFile(pars, template, outputpath=""):

	sys.stdout.write(pars["name"]+"...\n")
	formatpars = dict()
	for par in pars:
		formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	if (not os.path.exists(pars["outputpath"])):
		os.makedirs(pars["outputpath"])

	outfile = open(pars["outputpath"]+"/"+pars["name"]+".xml", 'w')
	outfile.write(output)
	outfile.close()
#


def formatPars(pars):

	# Parameter transformation
	if ("infectionRate" in pars.keys()):
	      #sys.stdout.write("Canonical parameterization detected\n\n")
	      
	      pars["becomeUninfectiousRate"] = pars["samplingRate"] + pars["recoveryRate"]
	      pars["reproductiveNumber"]     = float(pars["infectionRate"]) / pars["becomeUninfectiousRate"]
	      pars["samplingProportion"]     = float(pars["samplingRate"]) / pars["becomeUninfectiousRate"]

	elif ("reproductiveNumber" in pars.keys()):
	      #sys.stdout.write("reproductiveNumber parameterization detected\n\n")

	      pars["infectionRate"] = pars["reproductiveNumber"]*pars["becomeUninfectiousRate"]
	      pars["samplingRate"] 	= pars["samplingProportion"] * pars["becomeUninfectiousRate"]
	      pars["recoveryRate"]  = pars["becomeUninfectiousRate"] - pars["samplingRate"]            

	elif ("infectionRate1" in pars.keys() and "infectionRate2" in pars.keys()):
	      #sys.stdout.write("Canonical parameterization detected\n\n")
	      
	      pars["becomeUninfectiousRate"] = pars["samplingRate"] + pars["recoveryRate"]
	      pars["reproductiveNumber1"]    = float(pars["infectionRate1"]) / pars["becomeUninfectiousRate"]
	      pars["reproductiveNumber2"]    = float(pars["infectionRate2"]) / pars["becomeUninfectiousRate"]
	      pars["samplingProportion"]     = float(pars["samplingRate"]) / pars["becomeUninfectiousRate"]

	elif ("reproductiveNumber1" in pars.keys() and "reproductiveNumber2" in pars.keys()):
	      # sys.stdout.write("reproductiveNumber parameterization detected\n\n")

	      pars["infectionRate1"]  = pars["reproductiveNumber1"]*pars["becomeUninfectiousRate"]
	      pars["infectionRate2"]  = pars["reproductiveNumber2"]*pars["becomeUninfectiousRate"]
	      pars["samplingRate"] 	  = pars["samplingProportion"] * pars["becomeUninfectiousRate"]
	      pars["recoveryRate"]    = pars["becomeUninfectiousRate"] - pars["samplingRate"]            

	else :
	      sys.stderr.write("Unknown parameterization! exiting...\n")
	      sys.exit()
	#

	# Adjust for density dependence and remove first infected
	if ("popSize" in pars.keys()):
		pars["infectionRate"] = pars["infectionRate"]/pars["popSize"]
		pars["popSize"]   = pars["popSize"] - 1
	#
#


################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
	if (fnmatch(filename,config)):

		configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
		pars = yaml.load(configfile)

		# Set parameters not in the config file
		formatPars(pars)

		# Replace config file and save
		template = open(pars['template'], 'r').read()
		makeXMLFile(pars, template)
	#
#