#!/usr/bin/env python
'''
Created on Jun 8, 2011
A program to generate a dotgraph based on information acquired from a vcf file, using genezoom to create a table of said data
@author: jcc7
'''
import vcf
#from vcfutils import *
import re
import CrossTable
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Ellipse
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from optparse import OptionParser
import bed

######################################################

#Method to convert a base pair to its location as a base pair in the exon
#Returns -1 if the gene is not inside an exon
def bp2exonbp(tupleList, number):
    sum=0
    #whether it should be 0-based or 1-based
    base=0
    number=int(number)
    for entry in tupleList:
        if (number>=entry[0])&(number<entry[1]):
            return sum+number+1-entry[0]+base-1
        elif (number)<entry[0]:
            return -1
        else:
            sum=sum+entry[1]-entry[0]

#draw a exon at the proper location
def drawExon(exonTupleList):
    #Rectangle (xpos (left), ypos (bottom), width, height, kwargs)
    patches=[Rectangle((exonTupleList[0][0], -1), 1, 2, fill=False)]
    for exon in exonTupleList:
        start=exon[0]
        end=exon[1]
        patches.append(Rectangle((start, -1), end-start-1, 2))
    return PatchCollection(patches, match_original=True)

#draw stacked circles in a number equal to circles
def multiCircles(patches, circles, xLoc, circleLoc, circleWidth, circleHeight, colorShade):
    i=0
    #as long as there are circles left to draw, keep drawing them
    while i<circles:
        #Circle(size 2 array detailing xy, width)
        patches.append(Ellipse([xLoc, circleLoc], circleWidth, circleHeight, color=colorShade, linewidth=1.0))
        #if we are drawing downward, draw the next circle below this one, otherwise draw the next circle above
        if circleLoc<0:
            circleLoc=circleLoc-(circleHeight)
        else:
            circleLoc=circleLoc+(circleHeight)
        i+=1    
    return circleLoc #so that we know where to start drawing the next circles

#create a dotgraph, receiving a crosstable and an xLoc
def dotPlot(stuff, xLoc):
    #initialization of our patches, with a blank circle to avoid errors of empty list
    patches=[Circle([xLoc, 0], 1, color='white', alpha=0)] 
    circleWidth=0.4
    circleHeight=1.0
    for i in range(2):
        #Graph of the first list of data ('case') in our data set
        if i==0:
            circleLoc=circleHeight/2
        else:    
            circleLoc=-circleHeight/2
        #how many 0/1 alleles exist, and how many 1/1 alleles exist? This assumes that they are sorted onto the end of the column
        length=len(stuff.getbkeys())
        oneCircles= stuff.valueAt(i, length-2)
        twoCircles=stuff.valueAt(i, length-1)
        if twoCircles!=0:
            colorShade='#dd0000'
            circleLoc=multiCircles(patches, twoCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        if oneCircles!=0:
            colorShade='#00dd00'
            circleLoc=multiCircles(patches, oneCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
    return PatchCollection(patches, match_original=True)

def SetupPlot(start, end):
        #set up the graph format
        fig=plt.figure(facecolor='white')
        #add axes in rectangle left position, bottom position, width, height
        ax1 = fig.add_axes([0.1, 0.3, 0.8, 0.6])
        ax2=fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
        ax2.set_yticks([])
        ax1.set_title("Frequency of Alleles")
        ax1.set_ylabel("case                         control")
        ax2.set_xlabel("Chromosome")
        halfrange=(end-start)/2.0
        center=start+halfrange
        ax1.axis([center-halfrange, center+halfrange, -7, 7])
        horizontalLine=Line2D([-1, 100000000000], [0, 0],linewidth=1, color='black')
        ax1.add_line(horizontalLine)
        ax1.grid(True)
        ax2.set_ylim(-5, 5)
        ax2.add_line(horizontalLine)
        #fig.subplots_adjust(bottom=0.2)

        return ax1, ax2
    

############################################################
if __name__ == "__main__":
    #Begin option parser
    parser=OptionParser()
    parser.add_option(
        "-q", "--quiet",
        dest="verbose", 
        action="store_false", 
        default=True,
        help="don't print status messages to stdout"
        )

    parser.add_option(
        "-i", "--interact",
        dest="interact", 
        action="store_true", 
        default=False,
        help="Enter interactive python session after running"
        )

    parser.add_option(
        "-b", "--build",
        dest="build", 
        default='hg19',
        help="hg build (UCSC database for gene data)"
        )

    parser.add_option(
        "--bed",
        dest="bed", 
        default='../testing/data/refFlat.txt.gz.1',
        metavar="refFlat|knownGene|filename",
        help="UCSC table or file to use for gene information"
        )
    
    parser.add_option(
        "--gene",
        dest="gene",
        default='NM_152486',
        help="Gene to graph")

    parser.add_option(
        "-v", "--vcf", 
        dest="vcf_file", 
        default="../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1",
        help ="vcf file containing genotypes", 
        metavar="FILE"
        )

    parser.add_option(
        "-t", "--traits", 
        dest="trait_file", 
        default="../testing/data/458_traits.csv",
        help="trait file", 
        metavar="FILE"
        )

    parser.add_option(
        "-g", "--groups", 
        dest="groups", 
        default="T2D",
        help="specify grouping variable", 
        metavar="STRING"
        )

    parser.add_option(
        "-r", "--region", 
        dest="region", 
        default='1:873000-880000',
        help="specify region of interest", 
        metavar="chr:start-stop"
        )

    parser.add_option(
        "-p", "--prefix", 
        dest="prefix", 
        default="genezoom-out",
        help ="prefix for output files", 
        metavar="STRING"
        )

    (options, args) = parser.parse_args()

    #Begin regular expression compile for chrom, start, and stop
    regionRE=re.compile(r'(.+):(\d+)-(\d+)')
    
    m=regionRE.match(options.region)
    chrom=int(m.groups()[0])
    start=int(m.groups()[1])
    end=int(m.groups()[2])
    try:
        from vcfutils import *        
        v=vcfReader(options.vcf_file)
        vstuff = v.query(chrom, start, end)
    except Exception as e:
        print e
        v = vcf.VCFdata(options.vcf_file, 1024*1024)
        vstuff = v.fetch_range(chrom, start, end)
    
    refFlatKeys = ['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
    refFlat = bed.BED(options.bed, keys=refFlatKeys)
    bedRow=[row for row in refFlat if row['name']==options.gene ][0]
    
    #check for dimensions of graph, and to see if the introns intersect with the desired region
    ax1, ax2=SetupPlot(start, end)
    
    xLoc=start


    #for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
    
    for v in vstuff:
        #check to see if the gene is in the exon
        if not bp2exonbp(bedRow.get_exons(), v.getpos())==-1:
            xTable=CrossTable.xTable(options.trait_file, v.genotypes())
            dots=dotPlot(xTable, v.getpos())
            ax1.add_collection(dots)
    
    exonRect=drawExon(bedRow.get_exons())
    ax2.add_collection(exonRect)
    ax2.add_line(Line2D([-1, 100000000000], [0, 0],linewidth=1, color='black'))
    ax1.add_line(Line2D([-1, 100000000000], [0, 0],linewidth=1, color='black'))
    #ax1.set_xticks([xtickgathered])
    plt.show()

##################################################################################