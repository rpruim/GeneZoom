#!/usr/bin/env python

from gzutils import *
import gzio
from optparse import OptionParser, OptionGroup
import bed
import re
import genePlot as gp
import logging
import shlex
import sys

def parseCommandLine(line):
	my_splitter = shlex.shlex(line, posix=True)
	my_splitter.whitespace_split = True 
	return list(my_splitter)

# Load program constants.
conf_file = find_relative("conf/gz.conf")
if conf_file:
	execfile(conf_file)
else:
	die('Unable to find configuration file.')

def OptionSetUp(additional_args = ''):
	'''Sets up the option parser for the command line, returns the options chosen.'''
	#Begin option parser
	parser=OptionParser()
	#set up groups to consolidate the parser options
	graphGroup=OptionGroup(parser, "Graph options")
	infoGroup=OptionGroup(parser, "Genetic information options")
	outputGroup=OptionGroup(parser, "Output options")
	parser.add_option(
		"-b", "--batch",
		dest="batchfile", 
#		default=None,
		metavar="FILENAME",
		help="file with controls for batch operation")
	parser.add_option(
		"-q", "--quiet",
		dest="verbose", 
		action="store_false", 
		default=True,
		help="don't print status messages to stdout")
	parser.add_option(
		"-d", "--delim",
		dest="delim", 
		default=",",
		help="Delimiter for text files")
	parser.add_option(
		"-i", "--interact",
		dest="interact", 
		action="store_true", 
		default=False,
		help="Enter interactive python session after running")
	infoGroup.add_option(
		"", "--build",
		dest="build", 
		default='hg19',
		help="hg build (UCSC database for gene data)")
	infoGroup.add_option(
		"--gene",
		dest="gene",
		#default='NM_152486',
		metavar="GENENAME",
		help="Gene to graph")
	infoGroup.add_option(
		"-w", "--directory",
		dest = "directory",
		default = "./",
		help = "directory of data files")
	infoGroup.add_option(
		"--bed",
		dest="bed", 
		default='../testing/data/refFlat.txt.gz.1',
		metavar="refFlat|knownGene|filename",
		help="UCSC table or file to use for gene information")
	infoGroup.add_option(
		"-v", "--vcf", 
		dest="vcf_file", 
		default="../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1",
		help ="vcf file containing genotypes", 
		metavar="FILE")
	infoGroup.add_option(
		"-t", "--traits", 
		dest="trait_file", 
		default="../testing/data/458_traits.csv",
		help="trait file", 
		metavar="FILE")
	infoGroup.add_option(
		"-g", "--groups", 
		dest="groups", 
		default="T2D",
		help="specify grouping variable in trait file", 
		metavar="STRING")
	infoGroup.add_option(
		"--id",
		dest= "id",
		default = "ID",
		help = "specify ID variable in trait file")
	infoGroup.add_option(
		"--filter",
		dest = "filter",
		default = None,
		help = "Specify included quality filter in format filter1,filter2, etc.")
	graphGroup.add_option(
		"--title",
		dest="title",
		metavar='STRING',
		help="title for the plot")
	graphGroup.add_option(
		"-f", "--flank", 
		dest="flank", 
		default='0,0',
		help="flanking amount (bp)", 
		metavar="left[,right]")
	graphGroup.add_option(
		"-r", "--region", 
		dest="region", 
		#default='1:873000-880000',
		help="specify region of interest", 
		metavar="chr:start-stop")
	outputGroup.add_option(
		"--graph",
		dest="graph",
		default=True,
		action="store_true",
		help="Show the plot in an interactive session (default)")
	outputGroup.add_option(
		"--nograph",
		dest="graph",
		action="store_false",
		help="Don't show the plot in an interactive session")
	outputGroup.add_option(
		"-p", "--prefix", 
		dest="prefix", 
		default="results/genezoom",
		help ="filename for output files", 
		metavar="STRING")
	outputGroup.add_option(
		"--png",
		dest="png",
		action="store_true",
		help="Save graph to png")
	outputGroup.add_option(
		"--pdf",
		dest="pdf",
		action="store_true",
		help="Save graph to pdf")
	graphGroup.add_option(
		"-y", "--yscale", 
		dest="yscale", 
		default='25:25',
		help="number of controls shown: number of cases shown", 
		metavar="controls:cases")
	graphGroup.add_option(
		"--introns",
		dest="introns",
		action="store_true",
		help="Show the introns in the graph")
	graphGroup.add_option(
		"--nointrons",
		dest="introns",
		default=False,
		action="store_false",
		help="Don't show the introns in the graph (default)")
	graphGroup.add_option(
		"--codons",
		dest="codons",
		action="store_true",
		help="Graph data as codons.")
	graphGroup.add_option(
		"--nocodons",
		dest="codons",
		default=False,
		action="store_false",
		help="Graph data as nucleotides (default)")
	graphGroup.add_option(
		"--shape",
		dest="shape",
		default="circle",
		help="Plotted shapes.  Options: 'circle','rectangle' (default circle)")
	graphGroup.add_option(
		"-c", 
		"--color",
		dest="color",
		default="#0e51a7:#0acf00:#ff9e00:#fd0006",
		help = "RGB colors for graph.  List of colors in format: #012345:#6789AB:#2468AC:#13579B, for (1/1 frequency, 1/0 frequency, exonColor, exonColor)" )
	parser.add_option_group(infoGroup)
	parser.add_option_group(graphGroup)
	parser.add_option_group(outputGroup)
	(options, args) = parser.parse_args(sys.argv[1:] + parseCommandLine(additional_args))
	
	return (options, args)


def DataSetup( traitfile, bedfile ):
	'''Sets up the data for the UCSC file for gene information and the exon base pairs.'''
	traits = gzio.read_csv(traitfile, options.delim)
	refFlatKeys = ['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
	refFlat = bed.BED(bedfile, keys=refFlatKeys)
	return refFlat, traits

def PrintOptions( options, region ):
	'''Prints out the chosen options.  Can be used for debugging or user info.'''
	print "Chosen options:"
	print "  Genetic information options:"
	print "	hg build:\t\t%s"%options.build
	print "	Grouping:\t\t%s"%options.groups
	print "	Gene:\t\t\t%s"%options.gene 
	print " Directory of information files:\t%s"%options.directory
	print "	UCSC bed:\t\t%s"%options.bed
	print "	VCF gene file:\t%s"%options.vcf_file
	print " Filtering options: %s"%options.filter
	
	print "  Graph options:"
	print "	Title: \t\t%s"%options.title
	print "	Region:\t\t%s: "%options.chrom+"%s-"%options.start+"%s"%options.stop
	print " Calculated Region:\t%s: %s-%s"%(region[0], region[1], region[2])
	print " Flanks: %s"%options.flank
	print " Calculated Flanks: %s, %s"%(options.flankList[0], options.flankList[1])
	print "	Cases/Controls: \t%s"%options.ymin+"/%s"%options.ymax
	print "	Show Introns: \t%s"%options.introns
	print " Codons:\t%s"%options.codons
	print "	Colors: \t\t%s"%options.color
	
	print "  Output options:"
	print "	Show graph? \t%s"%options.graph
	print "	Output prefix: \t%s"%options.prefix
	print "	Save as png? \t%s"%options.png
	print "	Save as pdf? \t%s"%options.pdf

def ProcessBed(bedrow, introns, region):
	# bedRow=[row for row in refFlat if row['name']==options.gene ][0]
	'''Make an exon dictionary of the base pairs, defaulting to an exon if chosen start is in an exon.'''
	if bedrow==[]: #if a bedrow is empty (that is, there is no gene), then create a false exonDict
		length = region[2]-region[1]+1
		exonDict = dict((region[1]+i,i) for i in range(0, length))
	else:
		exonDict=gp.exonbplist(bedrow.get_exons(), introns)
#	if not exonDict.has_key(region[1]):
#		options.local_start=gp.bp2exonbp(bedrow.get_exons(), region[1], options.introns)
#	else:
#		options.local_start=exonDict[region[1]]
#	if not exonDict.has_key(region[2]):
#		options.local_stop=gp.bp2exonbp(bedrow.get_exons(), region[2], options.introns)
#	else:
#		options.local_stop=exonDict[region[2]]

	return bedrow, exonDict

def DetermineRegion( options, bedrow=None ):
	'''parse the region, returning default region if no region is given'''
	# parse options.flank
	if options.flank:
		options.flankList = [ int(x) for x in (options.flank).split(',') ]
		if len(options.flankList) == 1:
			(options.flankList).append(options.flankList[0])
		if len(options.flankList) != 2:
			options.flankList = [ 0, 0 ]
	if options.region:
		try:
			regionRE=re.compile(r'(.+):(\d+)-(\d+)')
			m=regionRE.match(options.region)
			options.chrom = ((m.groups()[0]).lstrip('chr'))
			options.start = int(m.groups()[1])
			options.stop = int(m.groups()[2])
			return( options.chrom, options.start, options.stop )
		except Exception as e:
			print >> sys.stderr, e
			die( 'Invalid region specification:  ' + str(options.region ))

	else:
		if not bedrow:
			logging.critical( 'Unable to make plot as no region was specified and no gene was located.' )
			return (None, None, None)
		options.chrom = bedrow['chrom']
		scaleRE=re.compile(r'(.+)')
		reg=scaleRE.match(options.chrom)
		options.chrom=reg.groups()[0].strip('chr')
		options.start = int(bedrow['txStart'])
		options.stop = int(bedrow['txEnd'])
	if options.flank:
		return( options.chrom, options.start-options.flankList[0], options.stop+options.flankList[1] )
	return (options.chrom, options.start, options.stop)
def parseChoices(options):
	#evaluate the scale choices
	try:
		scaleRE=re.compile(r'(\d+):(\d+)')
		y=scaleRE.match(options.yscale)
		options.ymin=int(y.groups()[0])
		options.ymax=int(y.groups()[1])
	except Exception as e:
		print >> sys.stderr, e
		print "Invalid yscale region.  Defaulting to 25 controls, 25 cases."
		options.ymin=25
		options.ymax=25
	#evaluate the color choices
	try:
		splitColor=re.split('#([A-Fa-f0-9]{6})', options.color)
		options.colorallele2='#'+splitColor[1]
		options.colorallele1='#'+splitColor[3]
		options.exoncolor1='#'+splitColor[5]
		options.exoncolor2='#'+splitColor[7]
	except Exception as e:
		print >> sys.stderr, e
		print "Invalid color scheme.  Defaulting to original colors."
		options.colorallele2='#0e51a7'
		options.colorallele1='#0acf00'
		options.exoncolor1='#ff9e00'
		options.exoncolor2='#fd0006'
	#Check for given title for graph.  If none, default to gene.
	if options.title==None:
		options.plotTitle=options.gene
		if options.gene==None:
			options.plotTitle=options.region
	else:
		options.plotTitle=options.title
	#Check for desired shape.  If none, default to circle.
	if ((options.shape!="circle") and (options.shape!="rect") and (options.shape!="rectangle")):
		print "Invalid shape %s chosen. Defaulting to circle."%options.shape
		options.shape="circle"
	#Check for  s
	if options.filter!=None:
		options.filterList = [f.strip() for f in (options.filter).split(',')]
	else:
		options.filterList = None
	return options
		
def RunJob(job_options, bedRows, v, traits, region):
	'''Load and run job'''
	#logging.critical(str(region))
	(bedrow, exonDict) = ProcessBed(bedRows, job_options.introns, region)
	vstuff = v.reg2vcf(region[0], int(region[1]), int(region[2]))
	vcfIDs = v.get_headers()[9:]
	print len(vstuff), "markers in region " + str(region[0]) + ":" + str(region[1]) + "-" + str(region[2])
	#PrintOptions(options, region)
	gp.pictograph(job_options, vstuff, exonDict, bedrow, traits, region, vcfIDs)
######################################################################################	
if __name__ == "__main__":
	options, args = OptionSetUp()
	#PrintOptions( options )
	#load the vcf file containing the genotypes
	#this is placed here instead of in DataSetup due to import restrictions

	if options.batchfile:
		batch = open(options.batchfile, 'r')
		jobs = batch.readlines()
		jobs = [j for j in jobs if j[0] != '#' and len(j) > 1]

	else:
		jobs = ['']

	last_vcf_file = None
	last_trait_file = None
	last_traitfile = None
	for job in jobs:
		(job_options, job_args) = OptionSetUp(job)
		#PrintOptions( job_options )
		if not job_options.vcf_file == last_vcf_file:
			try:
				from tabix import *
				v=tabixReader(job_options.directory + job_options.vcf_file)
				last_vcf_file = job_options.vcf_file
			except Exception as e:
				print e
				logging.critical (str(e)) 
				logging.critical("Unable to open vcf file: " + str(job_options.vcf_file) )
				logging.critical("\tSkipping.")
				continue
		print last_vcf_file
		if job_options.trait_file != last_trait_file:
			traitfile = job_options.directory + job_options.trait_file
			refFlat, traits = DataSetup(traitfile, job_options.bed)
			last_trait_file = job_options.trait_file
		bedRows = refFlat.get_rows(job_options.gene)
		print "Gene Flavors =", len(bedRows)
		if len(bedRows) < 1:
			job_options=parseChoices( job_options )
			region = DetermineRegion( job_options )
			RunJob(job_options, bedRows, v, traits, region)
		else:
			for bedrow in bedRows:
				print "\n", bedrow['name'], bedrow['geneName']
				job_options=parseChoices(job_options)
				region = DetermineRegion( job_options, bedrow )
				RunJob(job_options, bedrow, v, traits, region)
	if options.interact:
		try:
			from IPython.Shell import IPShellEmbed  # enter interactive ipython
			ipshell = IPShellEmbed([])
			ipshell()
		except:
			logging.critical("Unable to locate ipython, so shutting down without entering interactive mode.")
	else:
		exit(0)
