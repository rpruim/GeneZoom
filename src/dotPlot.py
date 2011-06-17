#!/usr/bin/env python
'''
Created on Jun 8, 2011
A program to generate a dotgraph based on information acquired from a vcf file, using genezoom to create a table of said data
@author: jcc7
'''
import vcf, re
import CrossTable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from optparse import OptionParser

######################################################

#Method to convert a base pair to its location as a base pair in the exon
def bp2exonbp(tupleList, number):
    sum=0
    for entry in tupleList:
        if (number>=entry[0])&(number<=entry[1]):
            return sum+number+1-entry[0]
        elif (number)<entry[0]:
            return -1
        else:
            sum=sum+entry[1]-entry[0]+1

#draw an ellipse of color color
#def drawCircle(patches, xLoc, yLoc, radius, colorShade):
#    c=Circle([xLoc, yLoc], radius, color=colorShade)
#    patches.append(c)
#    return patches

#draw stacked circles in a number equal to circles
def multiCircles(patches, circles, xLoc, circleLoc, circleRadius, colorShade):
    i=0
    #as long as there are circles left to draw, keep drawing them
    while i<circles:
        #Circle(size 2 array detailing xy, width)
        patches.append(Circle([xLoc, circleLoc], circleRadius, color=colorShade, linewidth=1.0))
        #patches=drawCircle(patches, xLoc, circleLoc, circleRadius, colorShade)
        #if we are drawing downward, draw the next circle below this one, otherwise draw the next circle above
        if circleLoc<0:
            circleLoc=circleLoc-(circleRadius*2)
        else:
            circleLoc=circleLoc+(circleRadius*2)
        i+=1    
    return circleLoc #so that we know where to start drawing the next circles

#create a dotgraph, receiving a crosstable and an xLoc
def dotPlot(stuff, xLoc):
    #initialization of our patches, with a blank circle to avoid errors of empty list
    patches=[Circle([xLoc, 0], 1, color='white', alpha=0)] 
    circleRadius=0.15
    for i in range(2):
        #Graph of the first list of data ('case') in our data set
        if i==0:
            circleLoc=circleRadius
        else:    
            circleLoc=-circleRadius
        #these two lines will probably have to be modified based on size of table, etc.
        #how many 0/1 alleles exist, and how many 1/1 alleles exist? This assumes that they are sorted onto the end of the column
        length=len(stuff.getbkeys())
        oneCircles= stuff.valueAt(i, length-2)
        twoCircles=stuff.valueAt(i, length-1)
        if twoCircles!=0:
            colorShade='#ff0000'
            circleLoc=multiCircles(patches, twoCircles, xLoc, circleLoc, circleRadius, colorShade)
        if oneCircles!=0:
            colorShade='#00ff00'
            circleLoc=multiCircles(patches, oneCircles, xLoc, circleLoc, circleRadius, colorShade)
    return PatchCollection(patches, match_original=True)

def SetupPlot(start, end):
        #set up the graph format
        fig=plt.figure()
        plt.title("Frequency of Alleles")
        plt.xlabel("Chromosome")
        plt.ylabel("case                         control")
        ax1 = fig.add_subplot(111)  
        halfrange=(end-start)/2.0
        center=start+halfrange
        #half of the total desired range shown

        plt.axis([center-halfrange, center+halfrange, -2, 2])
        horizontalLine=Line2D([0, 100000000000], [0, 0],linewidth=1, color='black')
        ax1.add_line(horizontalLine)
        ax1.grid(True)
        fig.subplots_adjust(bottom=0.2)

        return ax1
    

############################################################
if __name__ == "__main__":
    #Begin option parser
    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename", default='../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf', help ="data input file", metavar="input")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
    parser.add_option("-r", "--region", dest="region", default="1:157506300-157511350", help="region requested, chrom: start-stop format", metavar="region")
    
    (options, args)=parser.parse_args()
    #Begin regular expression compile for chrom, start, and stop
    regionRE=re.compile(r'(.+):(\d+)-(\d+)')
    filename=options.filename
    
    m=regionRE.match(options.region)
    v = vcf.VCFdata(filename, 1024*1024)
    chrom=int(m.groups()[0])
    start=int(m.groups()[1])
    end=int(m.groups()[2])
    
    vstuff = v.fetch_range(chrom, start, end)
    
    
    ax1=SetupPlot(start, end)

    #example trait information, should be replaced
    traitstuff=[]
    j=0
    for v in vstuff:
        if j%2==0:
            traitstuff.append('case')
        else:
            traitstuff.append('control')
        j+=1
    
    xLoc=start
    #for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
    for v in vstuff:
        xTable=CrossTable.xTable(traitstuff, v.genotypes())
        dots=dotPlot(xTable, v.getpos())
        ax1.add_collection(dots)
    plt.show()

##################################################################################