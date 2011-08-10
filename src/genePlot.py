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


def drawExon(exonTupleList, exonDict, exoncolors, introns):
	'''Draws an exon at the proper location.'''
	#Rectangle (xpos (left), ypos (bottom), width, height, kwargs)
	patches = [Rectangle((0, -1), 1, 2, fill=False)] #initialize patches with a blank rectangle to avoid errors
	loc = 0 
	exonColor = 1

	for exon in exonTupleList:
		#width assumes exclusive nature for exons: that is, exon (20,22) is [20,22) and a width 1 rectangle will be drawn
		if introns: #if drawing introns
			start = exonDict[exon[0]] #make the start the location given by the exonDict
			width = exonDict[exon[1]-1] - exonDict[exon[0]] #make the width equal to the width between the locations
		else: #if drawing only exons
			start = loc #make the start at our current location
			width = exon[1] - exon[0] - 1 #make the width our distance between pairs
		#for exon colors
		if exonColor % 2 == 1: col = exoncolors[0] #if odd, use exoncolor1
		else: col = exoncolors[1] #if even, use exoncolor2
			
		patches.append(Rectangle((start, -1), width, 2, color=col)) #draw the exon
		loc += width + 1 #move the current location
		exonColor += 1
	return PatchCollection(patches, match_original=True)


def multiPatch(patches, patchAmount, xLoc, patchLoc, colorShade, shape):
	'''Draw stacked shapes as needed (rectangle or circle based on user-specified option, triangle if location is an indel).'''
	if patchLoc <0: direction=-1
	else: direction=1
	#as long as there are patches left to draw, keep drawing them
	for _ in range(patchAmount):
		if (shape == "rect") or (shape == "rectangle"):
			patches.append(Rectangle((xLoc - 0.2, patchLoc - 0.5), 0.4, 1, color=colorShade))#Rectangle 
		elif shape=="triangle":
			patches.append(RegularPolygon([xLoc, patchLoc], 3, radius=0.67, orientation =180, color=colorShade))
		else: #if shape==circle
			patches.append(Ellipse([xLoc, patchLoc], 0.4, 1, color=colorShade, linewidth=1.0))#Circle(size 2 array detailing xy, width)
		patchLoc += direction#location for the next patch
	return patchLoc #return the location of the current group of patches, so we know where to draw the next group


def patchPlot(stuff, xLoc, colors, shape, keys):
	'''Create a graph of tabulation, receiving a crosstable and an xLoc.'''
	#initialization of our patches, with a blank circle to avoid errors of empty list
	patches = [Circle([xLoc, 0], 1, color='white', alpha=0)] 
	for ccEntry in stuff.keys():
		if ccEntry == keys[0]: #if we are drawing controls
			patchLoc = -0.5 #draw them in a downward direction
		else: #otherwise
			patchLoc = 0.5 #draw them in an upward direction
		oneCount = 0 #initialize our counts (required to avoid errors if our table does not have this data)
		twoCount = 0
		#get our counts from our crosstable
		if stuff[ccEntry].has_key('1/0'):
			oneCount = stuff[ccEntry]['1/0']
		if stuff[ccEntry].has_key('0/1'):
			oneCount = oneCount + stuff[ccEntry]['0/1']
		if stuff[ccEntry].has_key('1/1'):
			twoCount = stuff[ccEntry]['1/1']
		patchLoc = multiPatch(patches, twoCount, xLoc, patchLoc, colors[1], shape)#draw 1/1 alleles
		patchLoc = multiPatch(patches, oneCount, xLoc, patchLoc, colors[0], shape)#draw 0/1 alleles
	return PatchCollection(patches, match_original=True) #return our collection of patches


def SetupPlot(dimensions, alleleColor, title, chrom):
	'''Set up the parameters for the graph.  Receives dimensions [xmin, xmax, ymin, ymax], alleleColor [color1, color2], title for the graph, and chromosome number.'''
	fig = plt.figure(figsize=(8, 5))#set window size to width, height 
	#add axes in rectangle left position, bottom position, width, height
	ax1 = fig.add_axes([0.1, 0.3, 0.8, 0.6])
	ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
	ax1.axis(dimensions) #set up axis ranges
	ax2.set_ylim(-5, 5)
	#set up titles, labels, and ticks
	ax1.set_title(title)
	ax2.set_xlabel("Chromosome %s"%chrom)
	ax1.grid(True)		
	ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int((abs(x))))))#set ticks to absolute value
	ax2.set_yticks([]) #eliminate yticks from the 2nd axis
	for label in ax1.xaxis.get_ticklabels():
		label.set_fontsize(0)#eliminate the labels from our top plot without eliminating them from the bottom
	horizontalLine = Line2D([-10000000000, 10000000000], [0, 0], linewidth=1, color='black')#draw horizontal x-axis
	ax1.add_line(horizontalLine)
	leg1 = ax1.legend((Ellipse((0, 0), 1, 1, color=alleleColor[1]), Ellipse((0, 0), 1, 1, color=alleleColor[0])), ('1/1', '1/0 and 0/1'), shadow=True, fancybox=True, prop=font_manager.FontProperties(size=10))
	try:
		leg1.draggable(state=True, use_blit=True)#make the legend draggable
	except:
		pass

	return ax1, ax2, fig

def graphName(fileTitle, dirname, extension):
	'''Saves file under directory, first checking to see if the directory exists, then if the file already exists, creating a new directory or new name based on needs.'''
	if not os.path.isdir("./" + dirname + "/"): #if the directory doesn't exist, then create it
		os.mkdir("./" + dirname + "/")
	filename= fileTitle + extension
	i=1
	while os.path.isfile(filename): #if the filename already exists, then make a new name of the form filenamenumber, where the number has 4 digits
		filename=fileTitle + "%.4d"%i + extension
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
	
def pictograph(options, vstuff, exonDict, bedRow, traits, region, vcfIDs):
	'''Creates a plot based upon a set of options, vcf information, a list of exon tuples, a bed of UCSC genomes, and a list of traits.'''
	start, stop = tuplesDomain((region[1], region[2]), options.introns, exonDict, bedRow)#change the region specified to basepair values
	dimensions = (start, stop, options.ymin * -1, options.ymax)
	allelecolors=(options.colorallele1, options.colorallele2)
	exoncolors=(options.exoncolor1, options.exoncolor2)	
	ax1, ax2, fig = SetupPlot(dimensions, allelecolors, options.plotTitle, region[0]) #initialize the graph, with proper range and choices
	vstuffFiltered = [v for v in vstuff if v.checkFilter(options.filterList)]
	#for each element of vstuff (the data of chromosomes) create the cross table, add the proper dotGraph to the total plot
	tableKeys = ()
	for v in vstuffFiltered:
		#check to see if the gene is in the exon.  If it is, create a cross table, draw the dots and add them to the graph
		if (exonDict.has_key(int(v.get_pos()))):
			organizedList=CrossTable.cullList(vcfIDs, traits[options.id], traits[options.groups])#organize the traits into a list, returning a list of case/control/None corresponding to the vcfIDs
			xTable = CrossTable.xTable(organizedList, v.get_genotypes())
			if len(tableKeys)!=2: #check for keys in dictionary, for y-axis labels
				tableKeys = xTable.getTable().keys()
			if v.is_indel():  #if this gene is an indel, change shape to triangles
				drawings = patchPlot(xTable.getTable(), exonDict[int(v.get_pos())], allelecolors, "triangle", tableKeys)
			else:
				drawings = patchPlot(xTable.getTable(), exonDict[int(v.get_pos())], allelecolors, options.shape, tableKeys)
			ax1.add_collection(drawings)
			print xTable
	if tableKeys !=(): ax1.set_ylabel("%s                       %s"%(tableKeys[0], tableKeys[len(tableKeys)-1]))
	print "%s markers plotted after filtering."%len(vstuffFiltered)
	if bedRow!=[]:#as long as we're actually drawing exons (so a gene, not just a region)
		exonRect = drawExon(bedRow.get_exons(), exonDict, exoncolors, options.introns) #draw the exons, adding them to the plot
		ax2.add_collection(exonRect)
	ax2.add_line(Line2D([-10000000000, 10000000000], [0, 0], linewidth=1, color='black'))
	if options.png: #if user has chosen to save graph as a png, save it
		fig.savefig(graphName(options.prefix, "results", ".png"))
	if options.pdf: #if user has chosen to save graph as a pdf, save it
		fig.savefig(graphName(options.prefix, "results", ".pdf"))
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
