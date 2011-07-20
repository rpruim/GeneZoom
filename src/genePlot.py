#!/usr/bin/env python
'''
Created on Jun 8, 2011
A program to generate a dotgraph based on information acquired from a vcf file, using genezoom to create a table of said data

'''

import CrossTable
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Ellipse, RegularPolygon
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from matplotlib import font_manager

######################################################

def bp2exonbp(tupleList, number, introns):
    '''Method to convert a base pair to its location as a base pair in the exon.
    Receives a list of tuples and the number of the base pair, returns the location in the exons (or returns the last location in the previous exon).'''
    sum=0
    if introns:
        return number-tupleList[0][0]
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

def exonbplist(tupleList, introns):
    '''Receives the list of tuples, returns a dictionary with the exon bp locations accessed by keys pertaining to their base pair.'''
    exonDict={}
    location=0
    oldEntry=0    
    for entry in tupleList:
        if oldEntry !=0:
            location+=entry[0]-oldEntry-1
        for i in range(entry[0], entry[1]+1):
            exonDict[i]= location
            location +=1
        if introns:
            oldEntry=entry[1]
    return exonDict


def drawExon(exonTupleList, exonDict, options):
    '''Draws an exon at the proper location.'''
    #Rectangle (xpos (left), ypos (bottom), width, height, kwargs)
    patches=[Rectangle((0, -1), 1, 2, fill=False)]
    loc=0 
    exonColor=0

    for exon in exonTupleList:
        #width assumes inclusive nature for exons: that is, exon (20,22) is [20,22] and a width 2 rectangle will be drawn
        if options.introns:
            start=exonDict[exon[0]]
            width=exonDict[exon[1]]-exonDict[exon[0]]
        else:
            start=loc
            width=exon[1]-exon[0]
        if exonColor%2==1: col=options.exoncolor1
        else: col=options.exoncolor2 #for some extra color
            
        patches.append(Rectangle((start, -1), width, 2, color=col))
        loc=loc+width+1
        exonColor+=1
    return PatchCollection(patches, match_original=True)


def multiPatch(patches, patchAmount, xLoc, patchLoc, patchWidth, patchHeight, colorShade, shape):
    '''Draw stacked patchAmount as needed.'''
    i=0
    #as long as there are patches left to draw, keep drawing them
    while i<patchAmount:
        if (shape=="rect") or (shape=="rectangle"):
            patches.append(Rectangle((xLoc-patchWidth/2, patchLoc-patchHeight/2), patchWidth, patchHeight, color=colorShade))#Rectangle 
        else:
            patches.append(Ellipse([xLoc, patchLoc], patchWidth, patchHeight, color=colorShade, linewidth=1.0))#Circle(size 2 array detailing xy, width)
        #if we are drawing downward, draw the next shape below this one, otherwise draw the next shape above
        if patchLoc<0:
            patchLoc=patchLoc-(patchHeight)
        else:
            patchLoc=patchLoc+(patchHeight)
        i+=1    
    return patchLoc #so that we know where to start drawing the next patch


def patchPlot(stuff, xLoc, options):
    '''Create a histogram graph, receiving a crosstable and an xLoc.'''
    #initialization of our patches, with a blank circle to avoid errors of empty list
    patches=[Circle([xLoc, 0], 1, color='white', alpha=0)] 
    patchWidth=0.4
    patchHeight=1.0
    for ccEntry in stuff.keys():
        #Graph of the first list of data ('case') in our data set
        if ccEntry=='control':
            patchLoc=-patchHeight/2
        else:    
            patchLoc=patchHeight/2
        oneCount=0
        twoCount=0
        if stuff[ccEntry].has_key('1/0'):
            oneCount=stuff[ccEntry]['1/0']
        if stuff[ccEntry].has_key('0/1'):
            oneCount=oneCount+stuff[ccEntry]['0/1']
        if stuff[ccEntry].has_key('1/1'):
            twoCount=stuff[ccEntry]['1/1']
        if twoCount!=0:
            colorShade=options.colorallele2
            patchLoc=multiPatch(patches, twoCount, xLoc, patchLoc, patchWidth, patchHeight, colorShade, options.shape)
        if oneCount!=0:
            colorShade=options.colorallele1
            patchLoc=multiPatch(patches, oneCount, xLoc, patchLoc, patchWidth, patchHeight, colorShade, options.shape)
    return PatchCollection(patches, match_original=True)


def SetupPlot(start, end, ymin, ymax, options):
    '''Set up the parameters for the graph.'''
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
    ax1.set_ylabel("control                          case")
    ax2.set_xlabel("Chromosome")
    ax1.grid(True)        
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int((abs(x))))))#set ticks to absolute value
    ax2.set_yticks([])
    for label in ax1.xaxis.get_ticklabels():
        # label is a Text instance
        label.set_fontsize(0)
        #add horizontal axis and a legend
    horizontalLine=Line2D([-1, 100000000], [0, 0],linewidth=1, color='black')
    ax1.add_line(horizontalLine)
    leg1=ax1.legend((Ellipse((0,0), 1, 1, color='#0e51a7'), Ellipse((0,0), 1, 1, color='#0acf00')),('1/1', '1/0 and 0/1'), shadow=True, fancybox=True, prop=font_manager.FontProperties(size=10))
    leg1.draggable(state=True, use_blit=True)
    #fig.subplots_adjust(bottom=0.2)

    return ax1, ax2, fig

def checkDir():
    '''Checks to see if a results directory exists.  If it does not, it creates one.'''
    dirname="results"
    if not os.path.isdir("./" + dirname + "/"):
        os.mkdir("./" + dirname + "/")

def histogram(options, vstuff, exonDict, bedRow, traits):
    '''Creates a plot based upon a set of options, vcf information, a list of exon tuples, a bed of UCSC genomes, and a list of traits.'''
    ax1, ax2, fig=SetupPlot(options.start, options.stop, options.ymin*-1, options.ymax, options)
    #for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
    for v in vstuff:
        #check to see if the gene is in the exon.  If it is, create a cross table, draw the dots and add them to the graph
        if exonDict.has_key(int(v.get_pos()) ):
            xTable=CrossTable.xTable(traits['T2D'], v.get_genotypes())
            drawings=patchPlot(xTable.getTable(), exonDict[int(v.get_pos())], options)
            ax1.add_collection(drawings)            
    exonRect=drawExon(bedRow.get_exons(), exonDict, options) #draw the exons
    ax2.add_collection(exonRect)
    ax2.add_line(Line2D([-1, 100000000000], [0, 0],linewidth=1, color='black'))
    if options.png:
        checkDir()
        fig.savefig("results/"+options.prefix+'.png')
    if options.pdf:
        fig.savefig("results/"+options.prefix+'.pdf')
    if options.graph:
        plt.show()
############################################################
if __name__ == "__main__":
    tuples=((22, 24), (37, 40), (52, 53), (55, 60))
    print exonbplist(tuples, False)
    print exonbplist(tuples, True)
