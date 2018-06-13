import sys, os, time, string
from fnmatch import fnmatch
from optparse import OptionParser
import subprocess
from subprocess import Popen, PIPE



# Check log files from BEAST runs and summarize HPDs
################################################################################################################################
# Parameters
################################################################################################################################

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)


parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to store output files [required]")

parser.add_option("-B","--BEAST",
                  dest = "beast",
                  default = "",
                  metavar = "path",
                  help = "Beast 2 jar file [required]")

parser.add_option("-e","--ess",
                  dest = "ess",
                  default = "",
                  metavar = "path",
                  help = "ESS considered sufficient [required]")

parser.add_option("-b","--burnin",
                  dest = "burnin",
                  default = "",
                  metavar = "path",
                  help = "Percentage of log-files and tree-files to discard as burnin [required]")

parser.add_option("-p","--pattern",
                  dest = "pattern",
                  default = "*.xml",
                  metavar = "",
                  help = "Pattern to search for (in quotes) [default = %default]")

parser.add_option("-m","--minsamples",
                  dest = "minsamples",
                  default = "1000000",
                  metavar = "",
                  help = "Minimum number of samples [default = %default]")

(options,args) = parser.parse_args()

inputpath  = os.path.abspath(options.inputpath)+ "/"
outputpath = os.path.abspath(options.outputpath)+ "/"
beast      = os.path.abspath(options.beast)
pattern    = options.pattern
ess_sufficient = int(options.ess)
minsamples = int(options.minsamples)
b          = float(options.burnin)
if (b > 0 and b < 1):
	burnin = int(round(b*100))
else:
	burnin = int(round(b))

skipstats = ["*Times*"]


################################################################################################################################  

def analyseESS(inputfile, skipstats = [], ess_sufficient=200): 

	tol = 1e-10
	insufficient_ess = []

	statfile = open(inputfile, "r")
	statfile.readline()
	for line in statfile:
		parts = line.strip().split()
		stat_name = parts[0]
		stat_ess  = float(parts[8])
		stat_95upper = float(parts[6])
		stat_95lower = float(parts[5])

		skip = False
		for stat in skipstats:
			if (fnmatch(stat_name, stat)):
				skip = True

		if (not skip):
			#sys.stdout.write("\t%s\t%f\n" % (stat_name, stat_ess))
			if (stat_ess < ess_sufficient and stat_95upper - stat_95lower > tol):
				insufficient_ess.append("%s\t%f" % (stat_name, stat_ess))

	statfile.close()
	return(insufficient_ess)
#	


def runLogAnalyser(inputpath, filename, outfile, beast, burnin=10):

	start = time.time()

	here = os.getcwd()
	os.chdir(inputpath)

	loganalyser = "java -cp "+beast+" beast.util.LogAnalyser"
		
	loganalyser_cmd = "%s -b %d %s > %s" % (loganalyser, burnin, filename, outfile)


	print(loganalyser_cmd)
	#subprocess.check_call(loganalyser_cmd.split(" "))

	handle = Popen(loganalyser_cmd, stdout=None, stderr=PIPE, shell=True)

	err = handle.stderr.read()
	if (err != ""):
	    sys.stderr.write("\tWarning! Errors encountered!\n")
	    #sys.stderr.write(err)

	os.chdir(here)

	end = time.time()
	sys.stdout.write("\tLog file analysed ("+str(end-start)+" seconds)\n\n")
#



################################################################################################################################  
start = time.time()

if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

runs = dict()
stats = []
for xmlfile in sorted(os.listdir(inputpath), reverse=True):
	if (fnmatch(xmlfile,pattern)):
		filename = xmlfile[:xmlfile.rfind('.')]+"_127.log"
		print(filename)
		if (os.path.exists(inputpath+filename)):				
			logfile = open(inputpath+filename,'r')
			for line in logfile:
				pass
			samples = line[:line.find('\t')]

			runLogAnalyser(inputpath, filename, outputpath+filename+".stats", beast, burnin)

			insufficient_ess = analyseESS(outputpath+filename+".stats", skipstats, ess_sufficient)
			if (len(insufficient_ess) > 0):
				sys.stdout.write("\tInsufficient ESS values:\n")
				for s in insufficient_ess:
					sys.stdout.write("\t\t"+s+"\n")
			
			stats.append((filename[:filename.rfind('.')], int(samples), len(insufficient_ess)))
		else:
			stats.append((filename[:filename.rfind('.')], -1, -1))


#
statfile = open(outputpath+"ESS.stats",'w')
statfile.write("Logfile\tSamples\tInsufficient ESS\n")

convfile = open(outputpath+"Converged.stats",'w')
statfile.write("Logfile\tSamples\n")

for s in sorted(stats, key = lambda s : s[1]):
	if (s[1] < minsamples or s[2] > 0):
		statfile.write("\t".join(map(str,s))+"\n")
	
	if (s[2] == 0):
		convfile.write("\t".join(map(str,s[:2]))+"\n")
		
statfile.close()
convfile.close()
end = time.time()
print("Total time taken: "+str(end-start)+" seconds")