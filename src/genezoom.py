#!/usr/bin/env python

from gzutils import *
import gzio
from vcfutils import *
from optparse import OptionParser
import dotGraph
import logging

# Load program constants.
conf_file = find_relative("conf/gz.conf")
if conf_file:
	execfile(conf_file)
else:
	die('Unable to find configuration file.')


def SetUp():
	import re
	parser=OptionParser()
	parser.add_option(
		"-q", "--quiet",
		dest="verbose", 
		action="store_false", 
		default=True,
		help="don't print status messages to stdout"
		)

	parser.add_option(
		"-b", "--build",
		dest="build", 
		default='hg19',
		help="hg build (UCSC database for gene data)"
		)

	parser.add_option(
		"--genes",
		dest="genes", 
		default='refFlat',
		metavar="refFlat|knownGene|filename",
		help="UCSC table or file to use for gene information"
		)

	parser.add_option(
		"-i", "--interact",
		dest="interact", 
		action="store_true", 
		default=False,
		help="Enter interactive python session after running."
		)

	parser.add_option(
		"-v", "--vcf", 
		dest="vcf_file", 
		default="../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz",
		help ="vcf", metavar="FILE"
		)

	parser.add_option(
		"-t", "--traits", 
		dest="trait_file", 
		default="../testing/data/458_traits.csv",
		help="trait file", metavar="FILE"
		)

	parser.add_option(
		"-g", "--groups", 
		dest="groups", 
		default="T2D",
		help="grouping variable", metavar="STRING"
		)

	parser.add_option(
		"-r", "--region", 
		dest="region", 
		default='1:68000-70000',
		help="grouping variable", metavar="chr:start-stop"
		)

	(options, args) = parser.parse_args()

	try:
		#m = re.search( r'(.*):(\d+)-(\d+)', options.region)
		m = re.search( r'(.*):(.+)-(.+)', options.region)
		options.chrom = m.groups()[0]
		options.start = int(m.groups()[1])
		options.stop = int(m.groups()[2])
	except Exception as e:
		print >> sys.stderr, e
		die( 'Invalid region specification:  ' + options.region )

	return (options, args)

def main( options ):
	print >> sys.stderr, options
	traits = gzio.read_csv(options.trait_file)
	print >> sys.stderr, [a for a in traits]
	status = traits[options.groups]
	print >> sys.stderr, tally(status)
	vcf = tabix.Tabix( options.vcf_file )
	print >> sys.stderr, vcf
	results = vcf.query(options.chrom,options.start,options.stop)
	
	for row in results:
		g = xtally( status, extract_genotypes(row) )
		print >> sys.stderr, g
		print "=" * 60
	print extract_genotypes(row)[:10]

	import plotting
	fig = dotGraph.SetupPlot(options.start, options.stop)
	for row in results:
		crossTable=CrossTable(status, extract_genotypes(row) )
		dotGraph(crossTable, extract_pos(row))

	return { 'vcf': vcf, 'status': status, 'fig':fig }

if __name__ == "__main__":
	import tabix
	options, args = SetUp()
	session = main(options)
	if options.interact:
		try: 
			from IPython.Shell import IPShellEmbed  # enter interactive ipython
			ipshell = IPShellEmbed([])
			ipshell()
		except:
			logging.critical("Unable to locate ipython, so shutting down without entering interactive mode.")

  	else:
		exit(0)


