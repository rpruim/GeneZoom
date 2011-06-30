#!/usr/bin/env python
'''
Created on Jun 8, 2011
A program to generate a dotgraph based on information acquired from a vcf file, using genezoom to create a table of said data
@author: jcc7
'''

import CrossTable
import matplotlib.pylab as plt
from matplotlib.patches import Circle, Rectangle, Ellipse
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from matplotlib import font_manager

######################################################
#Method to convert a base pair to its location as a base pair in the exon
def bp2exonbp(tupleList, number):
    sum=0
    #whether it should be 0-based or 1-based
    base=0
    number=int(number)
    for entry in tupleList:
        if (number>=entry[0])&(number<entry[1]):
            return sum+number-entry[0]+base
        elif (number)<entry[0]:
            return sum-1
        else:
            sum=sum+entry[1]-entry[0]
    return sum +1

#receives the list of tuples, returns a dictionary with the exon bp locations accessed by keys pertaining to their base pair
def exonbplist(tupleList):
    exonDict={}
    location=0
    for entry in tupleList:
        for i in range(entry[0], entry[1]):
            exonDict[i]= location
            location +=1
    return exonDict

#draw a exon at the proper location
def drawExon(exonTupleList):
    #Rectangle (xpos (left), ypos (bottom), width, height, kwargs)
    patches=[Rectangle((0, -1), 1, 2, fill=False)]
    loc=0
    exonColor=0
    for exon in exonTupleList:
        start=loc
        width=exon[1]-exon[0]-1
        if exonColor%2==1: col='#ff9e00'
        else: col='#fd0006' #for some extra color
        patches.append(Rectangle((start, -1), width, 2, color=col))
        loc=loc+width+1
        exonColor+=1
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
            colorShade='#0e51a7'
            circleLoc=multiCircles(patches, twoCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        if oneCircles!=0:
            colorShade='#0acf00'
            circleLoc=multiCircles(patches, oneCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
    return PatchCollection(patches, match_original=True)

#set up the parameters for the graph
def SetupPlot(start, end, ymin, ymax, options):
        #set up the graph format
        fig=plt.figure(figsize=(8,5))#set window size to width, height 
        #add axes in rectangle left position, bottom position, width, height
        ax1 = fig.add_axes([0.1, 0.3, 0.8, 0.6])
        ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
        halfrange=(end-start)/2.0
        center=start+halfrange
        ax1.axis([center-halfrange, center+halfrange, ymin, ymax]) #set up axis ranges
        ax2.set_ylim(-5, 5)
        #set up titles, labels, and ticks
        ax1.set_title(options.title)
        ax1.set_ylabel("case                         control")
        ax2.set_xlabel("Chromosome")
        ax1.grid(True)        
        ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int((abs(x))))))#set ticks to absolute value
        ax2.set_yticks([])
        for label in ax1.xaxis.get_ticklabels():
            # label is a Text instance
            label.set_fontsize(0)
        #add horizontal axis and a legend
        horizontalLine=Line2D([-1, 100000000000], [0, 0],linewidth=1, color='black')
        ax1.add_line(horizontalLine)
        leg1=ax1.legend((Ellipse((0,0), 1, 1, color='#0e51a7'), Ellipse((0,0), 1, 1, color='#0acf00')),('1/1', '1/0 and 0/1'), shadow=True, fancybox=True, prop=font_manager.FontProperties(size=10))
        leg1.draggable(state=True, use_blit=True)
        #fig.subplots_adjust(bottom=0.2)

        return ax1, ax2, fig

def dotHistogram(options, vstuff, exonDict, bedRow):
    ax1, ax2, fig=SetupPlot(options.start, options.stop, options.ymin*-1, options.ymax, options)
    #for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
    for v in vstuff:
        #check to see if the gene is in the exon
        if exonDict.has_key(int(v.get_pos())):
            xTable=CrossTable.xTable(options.trait_file, v.get_genotypes())
            dots=dotPlot(xTable, exonDict[int(v.get_pos())])
            ax1.add_collection(dots)
    exonRect=drawExon(bedRow.get_exons()) #draw the exons
    ax2.add_collection(exonRect)
    ax2.add_line(Line2D([-1, 100000000000], [0, 0],linewidth=1, color='black'))
    if options.png:
        fig.savefig(options.prefix+'.png')
    if options.pdf:
        fig.savefig(options.prefix+'.pdf')
    if options.graph:
        plt.show()
############################################################
