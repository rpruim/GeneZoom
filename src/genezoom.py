#!/usr/bin/env python

from gzutils import *
import gzio
from optparse import OptionParser, OptionGroup
import bed
import re
import genePlot as gp
import logging

# Load program constants.
conf_file = find_relative("conf/gz.conf")
if conf_file:
    execfile(conf_file)
else:
    die('Unable to find configuration file.')

def OptionSetUp():
    '''Sets up the option parser for the command line, returns the options chosen.'''
    #Begin option parser
    parser=OptionParser()
    #set up groups to consolidate the parser options
    graphGroup=OptionGroup(parser, "Graph options")
    infoGroup=OptionGroup(parser, "Genetic information options")
    outputGroup=OptionGroup(parser, "Output options")
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
        "-b", "--build",
        dest="build", 
        default='hg19',
        help="hg build (UCSC database for gene data)")
    infoGroup.add_option(
        "--gene",
        dest="gene",
        default='NM_152486',
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
        help="specify grouping variable", 
        metavar="STRING")
    graphGroup.add_option(
        "--title",
        dest="title",
        default="",
        help="desired title for the plot")
    graphGroup.add_option(
        "-r", "--region", 
        dest="region", 
        default='1:873000-880000',
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
        help ="prefix for output files", 
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
        help="number of cases shown: number of controls shown", 
        metavar="cases:controls")
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
    (options, args) = parser.parse_args()
    
    #use regular expressions to evaluate the user's choices
    #evaluate the region choices
    try:
        regionRE=re.compile(r'(.+):(\d+)-(\d+)')
        m=regionRE.match(options.region)
        options.chrom = (m.groups()[0])
        options.start = int(m.groups()[1])
        options.stop = int(m.groups()[2])
    except Exception as e:
        print >> sys.stderr, e
        die( 'Invalid region specification:  ' + options.region )
    #evaluate the scale choices
    try:
        scaleRE=re.compile(r'(\d+):(\d+)')
        y=scaleRE.match(options.yscale)
        options.ymin=int(y.groups()[0])
        options.ymax=int(y.groups()[1])
    except Exception as e:
        print >> sys.stderr, e
        print "Invalid yscale region.  Defaulting to 25 cases, 25 controls."
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
    return (options, args)


def dataSetup( options ):
    '''Sets up the data for the UCSC file for gene information and the exon base pairs.'''
    traits = gzio.read_csv(options.trait_file)
    refFlatKeys = ['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
    refFlat = bed.BED(options.bed, keys=refFlatKeys)
    bedRow=[row for row in refFlat if row['name']==options.gene ][0]
    #make an exon dictionary of the base pairs, defaulting to an exon if chosen start is in an exon
    exonDict=gp.exonbplist(bedRow.get_exons(), options.introns)
    if not exonDict.has_key(options.start):
        options.start=gp.bp2exonbp(bedRow.get_exons(), options.start)
    else:
        options.start=exonDict[options.start]
    if not exonDict.has_key(options.stop):
        options.stop=gp.bp2exonbp(bedRow.get_exons(), options.stop)
    else:
        options.stop=exonDict[options.stop]
    return exonDict, bedRow, options, traits

def printOptions( options ):
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

if __name__ == "__main__":
    options, args = OptionSetUp()
    #load the vcf file containing the genotypes
    #this is placed here instead of in dataSetup due to import restrictions
    try:
        from vcfutils import *
        v=vcfReader(options.vcf_file)
        vstuff = v.query(options.chrom, options.start, options.stop)
    except Exception as e:
        import vcf
        print e
        v = vcf.VCFdata(options.vcf_file, 1024*1024)
        #vcfutils requires options.chrom to be a string, whereas vcf requires options.chrom to be an int
        vstuff = v.fetch_range(int(options.chrom), options.start, options.stop)
    exonDict, bedRow, options, traits = dataSetup(options)
    gp.histogram(options, vstuff, exonDict, bedRow, traits)
    if options.interact:
        try:
            from IPython.Shell import IPShellEmbed  # enter interactive ipython
            ipshell = IPShellEmbed([])
            ipshell()
        except:
            logging.critical("Unable to locate ipython, so shutting down without entering interactive mode.")
    else:
        exit(0)



