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
#        default=None,
        metavar="FILENAME",
        help="file with controls for batch operation")
    parser.add_option(
        "-q", "--quiet",
        dest="verbose", 
        action="store_false", 
        default=True,
        help="don't print status messages to stdout")
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
#        default='NM_152486',
        metavar="GENENAME",
        help="Gene to graph")
    infoGroup.add_option(
        "--bed",
        dest="bed", 
        default='../testing/data/refFlat.txt.gz.1',
        metavar="refFlat|knownGene|filename",
        help="UCSC table or file to use for gene information")
    infoGroup.add_option(
        "-v", "--vcf", 
        dest="vcf_file", 
#        default="../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1",
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
        help="specify grouping variable", 
        metavar="STRING")
    graphGroup.add_option(
        "--title",
        dest="title",
        default="",
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
#        default='1:873000-880000',
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
        default="genezoom",
        help ="filename for output files", 
        metavar="STRING")
    outputGroup.add_option(
        "--directory",
        dest = "directory",
        default = "results",
        help = "output directory")
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
        "--shape",
        dest="shape",
        default="circle",
        help="Plotted shapes.  Options: 'circle','rectangle' (default circle)")
    graphGroup.add_option(
        "-c", 
        "--color",
        dest="color",
        default="#0e51a7:#0acf00:#ff9e00:#fd0006",
        help = "RGB colors for graph.  List of colors in format: #012345:#6789AB:#2468AC:#13579B, for (1/1 frequency, 1/0 frequency, exonColor, exonColor)"                  )
    parser.add_option_group(infoGroup)
    parser.add_option_group(graphGroup)
    parser.add_option_group(outputGroup)
    (options, args) = parser.parse_args(sys.argv[1:] + parseCommandLine(additional_args))
    
    #use regular expressions to evaluate the user's choices
    #evaluate the region choices
#    try:
#        regionRE=re.compile(r'(.+):(\d+)-(\d+)')
#        m=regionRE.match(options.region)
#        options.chrom = (m.groups()[0])
#        options.start = int(m.groups()[1])
#        options.stop = int(m.groups()[2])
#    except Exception as e:
#        print >> sys.stderr, e
#        die( 'Invalid region specification:  ' + options.region )
#    #evaluate the scale choices
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
    if options.title=="":
        options.title=options.gene
    #Check for desired shape.  If none, default to circle.
    if ((options.shape!="circle") and (options.shape!="rect") and (options.shape!="rectangle")):
        print "Invalid shape %s chosen. Defaulting to circle."%options.shape        
        options.shape="circle"

	# parse options.flank
	if options.flank:
		options.flank = [ int(x) for x in (options.flank).split(',') ]
		if len(options.flank) == 1:
			(options.flank).append(options.flank[0])
		if len(options.flank) != 2:
			options.flank = [ 0, 0 ]

    return (options, args)


def DataSetup( options ):
    '''Sets up the data for the UCSC file for gene information and the exon base pairs.'''
    traits = gzio.read_csv(options.trait_file)
    refFlatKeys = ['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
    refFlat = bed.BED(options.bed, keys=refFlatKeys)
    return refFlat, traits

def PrintOptions( options ):
    '''Prints out the chosen options.  Can be used for debugging or user info.'''
    print "Chosen options:"
    print "  Genetic information options:"
    print "    hg build:\t\t%s"%options.build
    print "    Grouping:\t\t%s"%options.groups
    print "    Gene:\t\t%s"%options.gene 
    print "    UCSC bed:\t\t%s"%options.bed
    print "    VCF gene file:\t%s"%options.vcf_file
    
    print "  Graph options:"
    print "    Title: \t\t%s"%options.title
    print "    Region:\t\t%s: "%options.chrom+"%s-"%options.start+"%s"%options.stop
    print "    Cases/Controls: \t%s"%options.ymin+"/%s"%options.ymax
    print "    Show Introns: \t%s"%options.introns
    print "    Colors: \t\t%s"%options.color
    
    print "  Output options:"
    print "    Show graph? \t%s"%options.graph
    print "    Output prefix: \t%s"%options.prefix
    print "    Save as png? \t%s"%options.png
    print "    Save as pdf? \t%s"%options.pdf

def ProcessBed(bedrow, options):
    # bedRow=[row for row in refFlat if row['name']==options.gene ][0]

    #make an exon dictionary of the base pairs, defaulting to an exon if chosen start is in an exon
    exonDict=gp.exonbplist(bedrow.get_exons(), options.introns)
    if not exonDict.has_key(options.start):
        options.local_start=gp.bp2exonbp(bedrow.get_exons(), options.start, options.introns)
    else:
        options.local_start=exonDict[options.start]
    if not exonDict.has_key(options.stop):
        options.local_stop=gp.bp2exonbp(bedrow.get_exons(), options.stop, options.introns)
    else:
        options.local_stop=exonDict[options.stop]

    return bedrow, exonDict, options

def DetermineRegion( options, bedrow ):
    if options.region:
        try:
            regionRE=re.compile(r'(.+):(\d+)-(\d+)')
            m=regionRE.match(options.region)
            options.chrom = (m.groups()[0])
            options.start = int(m.groups()[1])
            options.stop = int(m.groups()[2])
            return( options.chrom, options.start, options.stop )
        except Exception as e:
            print >> sys.stderr, e
            die( 'Invalid region specification:  ' + str(options.region ))

    else:
        options.chrom = bedrow['chrom'] 
        options.start = bedrow['txStart'] 
        options.stop = bedrow['txEnd']
        return( options.chrom, options.start, options.stop )

if __name__ == "__main__":
    options, args = OptionSetUp()
    refFlat, traits = DataSetup(options)
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
    for job in jobs:
        (job_options, job_args) = OptionSetUp(job)
        if job_options.vcf_file != last_vcf_file:
            try:
                from tabix import *
                v=tabixReader(options.vcf_file)
                logging.debug('Using tabix.py')
            except Exception as e:
                print e
                die("Unable to open vcf file: " + str(options.vcf_file) )
        if options.trait_file != last_trait_file:
            refFlat, traits = DataSetup(options)

        for bedrow in refFlat.get_rows(options.gene):
            region = DetermineRegion( job_options, bedrow )
            (bedrow, exonDict, options) = ProcessBed(bedrow, job_options)
            vstuff = v.reg2vcf(job_options.chrom, int(job_options.start), int(job_options.stop))
            print len(vstuff), "markers in region " + job_options.chrom + ":" + str(job_options.start) + "-" + str(job_options.stop)
            gp.pictograph(job_options, vstuff, exonDict, bedrow, traits)
        last_vcf_file = job_options.vcf_file
        last_trait_file = job_options.trait_file

    if options.interact:
        try:
            from IPython.Shell import IPShellEmbed  # enter interactive ipython
            ipshell = IPShellEmbed([])
            ipshell()
        except:
            logging.critical("Unable to locate ipython, so shutting down without entering interactive mode.")
    else:
        exit(0)



