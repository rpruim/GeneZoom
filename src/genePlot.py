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
    sum = 0
    if introns: #if we are graphing introns, then return the distance from our number to the beginning intron
        return number - tupleList[0][0]
    base = 0 #whether it should be 0-based or 1-based
    number = int(number)
    #for each tuple, check to see if our number is contained within the pair. If
    #it is, add our sum to where the number is in the pair.  If our number
    #existed before the current pair, then return one less than the front of the
    #pair.  Otherwise, add the amount in the pair and move on.  If we reach the
    #end of the tupleList, then place our location at the end of the tupleList.
    for entry in tupleList:
        if (number >= entry[0]) & (number < entry[1]):
            return sum + number - entry[0] + base
        elif (number) < entry[0]:
            return sum - 1
        else:
            sum = sum + entry[1] - entry[0]
    return sum-1

def exonbplist(tupleList, introns):
    '''Receives the list of tuples, returns a dictionary with the exon bp locations accessed by keys pertaining to their base pair.  Each element of the dictionary will be in the form basepair:location'''
    exonDict = {}
    location = 0
    oldEntry = 0    
    for entry in tupleList:
        if oldEntry != 0: #if we are not at the beginning of the list, and are graphing introns
            location += entry[0] - oldEntry #add the distance between the last entry and this one to the location
        for i in range(entry[0], entry[1]): #for each number within our tuple
            exonDict[i] = location #set its location to where it is in the tuple
            location += 1 #advance to next location
        if introns: #if we are graphing introns
            oldEntry = entry[1] #set oldEntry for location calculation
    return exonDict


def drawExon(exonTupleList, exonDict, options):
    '''Draws an exon at the proper location.'''
    #Rectangle (xpos (left), ypos (bottom), width, height, kwargs)
    patches = [Rectangle((0, -1), 1, 2, fill=False)] #initialize patches with a blank rectangle to avoid errors
    loc = 0 
    exonColor = 1

    for exon in exonTupleList:
        #width assumes exclusive nature for exons: that is, exon (20,22) is [20,22) and a width 1 rectangle will be drawn
        if options.introns: #if drawing introns
            start = exonDict[exon[0]] #make the start the location given by the exonDict
            width = exonDict[exon[1]-1] - exonDict[exon[0]] #make the width equal to the width between the locations
        else: #if drawing only exons
            start = loc #make the start at our current location
            width = exon[1] - exon[0] - 1 #make the width our distance between pairs
        #for exon colors
        if exonColor % 2 == 1: col = options.exoncolor1 #if odd, use exoncolor1
        else: col = options.exoncolor2 #if even, use exoncolor2
            
        patches.append(Rectangle((start, -1), width, 2, color=col)) #draw the exon
        loc = loc + width + 1 #move the current location
        exonColor += 1
    return PatchCollection(patches, match_original=True)


def multiPatch(patches, patchAmount, xLoc, patchLoc, patchWidth, patchHeight, colorShade, shape):
    '''Draw stacked shapes as needed.'''
    i = 0
    #as long as there are patches left to draw, keep drawing them
    while i < patchAmount:
        if (shape == "rect") or (shape == "rectangle"):
            patches.append(Rectangle((xLoc - patchWidth / 2, patchLoc - patchHeight / 2), patchWidth, patchHeight, color=colorShade))#Rectangle 
        else: #if shape==circle
            patches.append(Ellipse([xLoc, patchLoc], patchWidth, patchHeight, color=colorShade, linewidth=1.0))#Circle(size 2 array detailing xy, width)
        
        if patchLoc < 0:#if we are drawing downward,
            patchLoc = patchLoc - (patchHeight)#draw the next shape below the current shape
        else: #otherwise
            patchLoc = patchLoc + (patchHeight) #draw the next shape above the current shape
        i += 1#advance the patch counter
    return patchLoc #return the location of the current group of patches, so we know where to draw the next group


def patchPlot(stuff, xLoc, options):
    '''Create a graph of tabulation, receiving a crosstable and an xLoc.'''
    #initialization of our patches, with a blank circle to avoid errors of empty list
    patches = [Circle([xLoc, 0], 1, color='white', alpha=0)] 
    patchWidth = 0.4
    patchHeight = 1.0
    for ccEntry in stuff.keys():
        if ccEntry == 'control': #if we are drawing controls
            patchLoc = -patchHeight / 2 #draw them in a downward direction
        else: #otherwise
            patchLoc = patchHeight / 2 #draw them in an upward direction
        oneCount = 0 #initialize our counts (required to avoid errors if our table does not have this data)
        twoCount = 0
        #get our counts from our crosstable
        if stuff[ccEntry].has_key('1/0'):
            oneCount = stuff[ccEntry]['1/0']
        if stuff[ccEntry].has_key('0/1'):
            oneCount = oneCount + stuff[ccEntry]['0/1']
        if stuff[ccEntry].has_key('1/1'):
            twoCount = stuff[ccEntry]['1/1']
        #if there are 1/1 alleles to draw, then call multiPatch to draw them
        if twoCount != 0:
            colorShade = options.colorallele2
            patchLoc = multiPatch(patches, twoCount, xLoc, patchLoc, patchWidth, patchHeight, colorShade, options.shape)
        #if there are 0/1 alleles to draw, then call multiPatch to draw them
        if oneCount != 0:
            colorShade = options.colorallele1
            patchLoc = multiPatch(patches, oneCount, xLoc, patchLoc, patchWidth, patchHeight, colorShade, options.shape)
    return PatchCollection(patches, match_original=True) #return our collection of patches


def SetupPlot(dimensions, alleleColor, title):
    '''Set up the parameters for the graph.'''
    #set up the graph format
    fig = plt.figure(figsize=(8, 5))#set window size to width, height 
    #add axes in rectangle left position, bottom position, width, height
    ax1 = fig.add_axes([0.1, 0.3, 0.8, 0.6])
    ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
    ax1.axis(dimensions) #set up axis ranges
    ax2.set_ylim(-5, 5)
    #set up titles, labels, and ticks
    ax1.set_title(title)
    ax1.set_ylabel("control                          case")
    ax2.set_xlabel("Chromosome")
    ax1.grid(True)        
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int((abs(x))))))#set ticks to absolute value
    ax2.set_yticks([]) #eliminate yticks from the 2nd axis
    for label in ax1.xaxis.get_ticklabels():
        # label is a Text instance
        label.set_fontsize(0)
        #add horizontal axis and a legend
    horizontalLine = Line2D([-10000000000, 10000000000], [0, 0], linewidth=1, color='black')
    ax1.add_line(horizontalLine)
    leg1 = ax1.legend((Ellipse((0, 0), 1, 1, color=alleleColor[1]), Ellipse((0, 0), 1, 1, color=alleleColor[0])), ('1/1', '1/0 and 0/1'), shadow=True, fancybox=True, prop=font_manager.FontProperties(size=10))
    try:
        leg1.draggable(state=True, use_blit=True)#make the legend draggable
    except:
        pass

    return ax1, ax2, fig

def saveGraph(fileTitle, dirname, extension):
    '''Saves file under directory, first checking to see if the directory exists, then if the file already exists, creating a new directory or new name based on needs.'''
    if not os.path.isdir("./" + dirname + "/"): #if the directory doesn't exist, then create it
        os.mkdir("./" + dirname + "/")
    filename= dirname + "/" + fileTitle + extension
    i=1
    while os.path.isfile(filename): #if the filename already exists, then make a new name of the form filename(number)
        filename=dirname + "/" + fileTitle + "(%s)"%i + extension
        i+=1
    return filename #save the figure under our created filename
    
def tuplesDomain((opStart, opStop), introns, exonDict, bedRow):
    '''Finds the required x scale based upon the base pair location of the region given by the user.'''
    if not exonDict.has_key(opStart):
        start=bp2exonbp(bedRow.get_exons(), opStart, introns)
    else:
        start=exonDict[opStart]
    if not exonDict.has_key(opStop):
        stop=bp2exonbp(bedRow.get_exons(), opStop, introns)
    else:
        stop=exonDict[opStop]
    return start, stop
    
def pictograph(options, vstuff, exonDict, bedRow, traits, region):
    '''Creates a plot based upon a set of options, vcf information, a list of exon tuples, a bed of UCSC genomes, and a list of traits.'''
    start, stop = tuplesDomain((region[1], region[2]), options.introns, exonDict, bedRow)#change the region specified to basepair values
    dimensions = (start, stop, options.ymin * -1, options.ymax)
    colors=(options.colorallele1, options.colorallele2)
    ax1, ax2, fig = SetupPlot(dimensions, colors, options.title) #initialize the graph, with proper range and choices
    #for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
    for v in vstuff:
        #check to see if the gene is in the exon.  If it is, create a cross table, draw the dots and add them to the graph
        if exonDict.has_key(int(v.get_pos())):
            xTable = CrossTable.xTable(traits['T2D'], v.get_genotypes())
            drawings = patchPlot(xTable.getTable(), exonDict[int(v.get_pos())], options)
            ax1.add_collection(drawings)            
    exonRect = drawExon(bedRow.get_exons(), exonDict, options) #draw the exons
    ax2.add_collection(exonRect) #add the collections to the graph
    ax2.add_line(Line2D([-10000000000, 10000000000], [0, 0], linewidth=1, color='black'))
    if options.png: #if user has chosen to save graph as a png, save it
        fig.savefig(saveGraph(options.prefix, options.directory, ".png"))
    if options.pdf: #if user has chosen to save graph as a pdf, save it
        fig.savefig(saveGraph(options.prefix, options.directory, ".pdf"))
    if options.graph: #if user has chosen to show the graph, then show it
        plt.show()
############################################################
if __name__ == "__main__":
    tupleList=((3, 6), (10, 14), (16, 19))
    print tupleList
    print "{",
    for i in range(3, 20):
        print "%s: %s,"%(i, bp2exonbp(tupleList, i, False)),
    print "}"
    print exonbplist(tupleList, False)
    print exonbplist(tupleList, True)
