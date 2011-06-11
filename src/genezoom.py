#!/usr/bin/env python

import gzio
from gzutils import *
from vcfutils import *
from optparse import OptionParser

# Load program constants.
conf_file = find_relative("conf/gz.conf");
execfile(conf_file);


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

if __name__ == "__main__":
	import tabix
	options, args = SetUp()
	main(options)

