#!/usr/bin/env python
'''
Created on Jun 8, 2011
A program to generate a dotgraph based on information acquired from a vcf file, using genezoom to create a table of said data
@author: jcc7
'''
from vcfutils import *
import sys
import tabix
import re
import CrossTable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from optparse import OptionParser

######################################################
#draw an ellipse of color color
def drawEllipse(ax1, xLoc, yLoc, width, height, color):
    e=Ellipse([xLoc, yLoc], width, height)
    e.set_facecolor(color)
    ax1.add_artist(e)

#draw stacked ellipses in a number equal to circles
def multiEllipses(ax1, circles, xLoc, circleLoc, circleWidth, circleHeight, colorShade):
    i=0
    #as long as there are circles left to draw, keep drawing them
    while i<circles:
        #Ellipse(size 2 array detailing xy, width, height)
        drawEllipse(ax1, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        #if we are drawing downward, draw the next circle below this one, otherwise draw the next circle above
        if circleLoc<0:
            circleLoc=circleLoc-circleHeight
        else:
            circleLoc=circleLoc+circleHeight
        i+=1    
    return circleLoc #so that we know where to start drawing the next circles

#create a dotgraph, receiving a crosstable and an xLoc
def dotGraph(ax1, stuff, xLoc):
    #size is the scale of the ellipse/circle
    size=1
    #circleWidth and circleHeight are the basic sizes of the ellipses (not yet properly scaled)
    circleWidth=1*size
    circleHeight=.2*size
    for i in range(2):
        #Graph of the first list of data ('case') in our data set
        if i==0:
            circleLoc=circleHeight/2
        else:    
            circleLoc=-circleHeight/2
        #these two lines will probably have to be modified based on size of table, etc.
        #how many 0/1 alleles exist, and how many 1/1 alleles exist? This assumes that they are sorted onto the end of the column
        length=len(stuff.getbkeys())
        oneCircles= stuff.valueAt(i, length-2)
        twoCircles=stuff.valueAt(i, length-1)
        if twoCircles!=0:
            colorShade=[0,0,1]
            circleLoc=multiEllipses(ax1, twoCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        if oneCircles!=0:
            colorShade=[0,1,0]
            circleLoc=multiEllipses(ax1, oneCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)

def SetupPlot(start, end):
    #set up the graph format
    fig=plt.figure()
    plt.title("Frequency of Alleles")
    plt.xlabel("Chromosome")
    ax1 = fig.add_subplot(111)
    halfrange=(end-start)/2.0
    center=start+ halfrange
    #half of the total desired range shown
    plt.axis([center-halfrange, center+halfrange, -2, 2])
    ax1.add_line(Line2D([0, 100000000000], [0, 0],linewidth=1, color='black'))
    ax1.grid(True)
    fig.subplots_adjust(bottom=0.2)

    #xlabels=[]
    #for c in vstuff:
        #xlabels.append(c.getpos())
    #xrange=np.arange(center-30, center+30, 1)
    #plt.xticks(xrange, xlabels, rotation=65)
    y=[-1,0,1]
    plt.yticks(y,('case', '','control'), rotation=90)
    return ax1

#############################################################
if __name__ == "__main__":
    #Begin option parser
    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename", 
			default='../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz', 
			help ="data input file", metavar="FILE")

    parser.add_option("-q", "--quiet", 
			action="store_false", 
			dest="verbose", 
			default=True, 
			help="don't print status messages to stdout")

    parser.add_option("-r", "--region", dest="region", 
			default="1:157508000-157511000", 
			help="region requested, chrom: start-stop format", metavar="REGION")
    
    (options, args)=parser.parse_args()

    #Begin regular expression compile for chrom, start, and stop
    regionRE=re.compile(r'(.+):(\d+)-(\d+)')
    filename=options.filename

    m=regionRE.match(options.region)
    chrom=str(m.groups()[0])
    start=int(m.groups()[1])
    end=int(m.groups()[2])

    vcf = tabix.Tabix(filename)
    results = [row for row in vcf.query(chrom, start, end)]
   
    ax1 = SetupPlot(start, end)
    
    #insert method to get vector from traits file here
    #generate some example data, as reading traits isn't yet working for me.  Replace this with real trait info.
    traitstuff=[]
    j=0
    for row in results:
        if j%2==0:
            traitstuff.append('case')
        else:
            traitstuff.append('control')
        j+=1
    
    xLoc=start
    #for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
    print >> sys.stderr, len(results)
    for row in results:
        xTable=CrossTable.CrossTable(traitstuff, extract_genotypes(row))
        dotGraph(ax1, xTable, extract_pos(row))
        print >> sys.stderr, extract_pos(row)
    plt.show()
	
##################################################################################
