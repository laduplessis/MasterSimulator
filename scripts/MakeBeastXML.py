import os, sys, io, yaml, json
import numpy as np
from fnmatch import fnmatch
from optparse import OptionParser
from Bio import Phylo, Nexus

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

parser.add_option("-x","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Path to template XML file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-n","--name",
                  dest = "name",
                  default = "",
                  metavar = "path",
                  help = "Name of the runs [required]")

parser.add_option("-s","--seed",
                  dest = "seed",
                  default = "127",
                  metavar = "integer",
                  help = "Seed or comma separated list of seeds to use [required]")

parser.add_option("-q","--queue",
                  dest = "queue",
                  default = "4",
                  metavar = "integer",
                  help = "Queue to submit to [required]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"

################################################################################################################################  

def writeDateTrait(dates,  output):

	traits = []	
	for i in range(0,len(dates)):
		traits.append('\n\t\t\t\t\t%d=%.10f' %(i+1, dates[i]))
	output.write(','.join(traits)+'\n')
#

def writeAlignment(n, output):

	for i in range(0,int(n)):
		output.write('\n\t\t\t\t\t<sequence spec="Sequence" taxon="%d" value="?" />' % (i+1))
#

def writeReactions(model, output):

	for reaction in model['reactions']:
		output.write('\n\t\t\t\t\t %s' % reaction)

#

# Read sampling times from tree and MASTER simulation and check that they are the same before returning
# 
# If forwards = True  gives times as time from tMRCA
# If forwards = False gives times as time from most recent sample
#
# If fromOrigin and forwards = True then times are measured from the origin of the simulation
def getSamplingTimes(simtrajectories, newicktree, tol=1e-6, forwards=False, fromOrigin=False, offset=None):

	#print(simtrajectories['Y'])
	
	# Sampling times from json file
	nsamples      = int(simtrajectories['Y'][-1])
	samplingtimes = np.zeros(nsamples)
	sampling = np.diff(np.array(simtrajectories['Y']))
	j = len(samplingtimes)-1
	for i in range(0,len(sampling)):
		if (sampling[i] == 1):
			samplingtimes[j] = simtrajectories['t'][i+1]
			j -= 1
	#	

	# Origin is the time of the most recent sample, not the simulation end time! 
	# (MASTER simulations can continue after the last sample has been taken)
	# origin = max(simtrajectories['t'])
	origin = samplingtimes[0]
	samplingtimes = (origin-samplingtimes)
	
	# Sampling times from tree
	tree    = Nexus.Trees.Tree(newicktree)
	leaves  = tree.get_terminals()
	labels  = []
	heights = np.zeros(len(leaves))	
	for i in range(0,len(leaves)):
		label = int(tree.get_taxa(node_id=leaves[i])[0])
		heights[label-1] = tree.sum_branchlength(node=leaves[i])
	#
	tmrca = max(heights)
	treetimes = tmrca-heights


	#print(forwards)
	#print("%.3f\t%.3f\t%.3f\n" % (samplingtimes[0], tmrca, origin))
	#print(sorted(samplingtimes))
	#print(sorted(treetimes))
	
	# Compare times	
	i = 0
	for t in sorted(treetimes):
		#sys.stdout.write('\n\t\t\t%.13f=%.13f' %(t, samplingtimes[i]))
		if (samplingtimes[i]-t > tol):
			sys.stdout.write("ERROR! Sampling times do not agree! (%.5f vs %.5f)\n" % (samplingtimes[i], t))
			#sys.stdout.write("%f\t%f\n" % (samplingtimes[i], t))
			sys.exit()
		i += 1
	#


	if (forwards == True):
		treetimes = (tmrca-treetimes)

		# Gap between origin and tMRCA 
		if (fromOrigin == True):
			treetimes = treetimes + (origin-tmrca)


	# Process date offset
	if (offset != None):
		treetimes = offset + treetimes



	return(treetimes)
#

def makeXMLFile(pars, template, outputpath=""):

	sys.stdout.write(pars["name"]+"...\n")
	formatpars = dict()
	for par in pars:
		formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	if (outputpath == ""):
		outputpath = pars['outputpath']

	if (not os.path.exists(outputpath)):
		os.makedirs(outputpath)

	outfile = open(outputpath+"/"+pars["name"]+".xml", 'w')
	outfile.write(output)
	outfile.close()
#


def formatPars(pars, tree, model, simtrajectories):

	# Set inputtree, sampling times and dummy alignment
	output_align = io.StringIO()
	output_dates = io.StringIO()
	output_model = io.StringIO()

	# Add default date direction if not defined
	if ("dateTrait" not in pars.keys()):
		pars["dateTrait"] = "date"

	# Add default fromOrigin if not defined
	if ("fromOrigin" not in pars.keys()):
		pars["fromOrigin"] = False

	# Add default dateOffset if not defined
	if ("dateOffset" not in pars.keys()):
		pars["dateOffset"] = None

	samplingTimes = getSamplingTimes(simtrajectories, tree, forwards=(pars["dateTrait"] != "date-backward"), 
															fromOrigin=(pars["fromOrigin"]), offset=pars["dateOffset"])

	

	writeDateTrait(samplingTimes, output_dates)
	writeAlignment(len(samplingTimes), output_align)
	writeReactions(model, output_model)

	pars["tree"]   	  = tree
	pars["dates"] = output_dates.getvalue()
	pars["alignment"] = output_align.getvalue()
	pars["reactions"] = output_model.getvalue()
	pars["origin"]    = max(simtrajectories['t'])

	if ('incubationRate' in pars.keys() and 'becomeUninfectiousRate' in pars.keys()):
		pars['becomeUninfectiousRate'] = pars['incubationRate']*pars['becomeUninfectiousRate']/(pars['incubationRate']+pars['becomeUninfectiousRate'])
#


################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
	if (fnmatch(filename,config)):

		sys.stdout.write(filename+"\t"+config+"\n")

		# Load config file
		configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
		pars 	   = yaml.load(configfile)

		# Load MASTER simulation results
		treefile   = open(pars["outputpath"]+pars["name"]+".newick", 'r')
		jsonfile   = open(pars["outputpath"]+pars["name"]+".json", 'r')
		simdata    = json.load(jsonfile)		

		# Set BEAST specific parameters
		queue      = int(options.queue)
		seeds      = list(map(int, options.seed.split(',')))
		basename   = pars["name"] if options.name == '' else pars["name"]+"_"+options.name
		outputpath = os.path.abspath(pars["outputpath"] if options.outputpath == '' else options.outputpath)
		template   = open(os.path.abspath(pars["template"] if options.template == '' else options.template), 'r').read()

		# Output scripts
		if (not os.path.exists(outputpath)):
			os.makedirs(outputpath)
		scriptfile = open(outputpath+"/"+basename+".sh",'w')
		eulerfile  = open(outputpath+"/"+basename+".euler.sh",'w')

		i = 0
		for tree in treefile:
				
			# Set parameters not in the config file
			formatPars(pars, tree, simdata['sim']['model'], simdata['trajectories'][i])

			# Replace config file and save
			pars["name"] = basename+"_"+str(i)
			makeXMLFile(pars, template, outputpath=outputpath+'/'+basename)

			# Write command to scripts
			for seed in seeds:
				cmd ="java -jar -Xms2G -Xmx4G $1 -overwrite -seed %d %s" % (seed, pars["name"]+".xml") 
				scriptfile.write("cd %s\n%s\ncd ..\n" % (basename, cmd))
				eulerfile.write("cd %s\nbsub -W%s:0 -n 1a -o %s.euler.out -R 'rusage[mem=4096]' %s\ncd ..\n" % (basename, queue, pars["name"], cmd))

			i += 1
		#

		treefile.close()
		jsonfile.close()
		scriptfile.close()
		eulerfile.close()
	#
#
