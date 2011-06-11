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
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from optparse import OptionParser

######################################################
#draw an ellipse of color color
def drawEllipse(xLoc, yLoc, width, height, color):
    e=Ellipse([xLoc, yLoc], width, height)
    e.set_facecolor(color)
    ax1.add_artist(e)
#draw stacked ellipses in a number equal to circles
def multiEllipses(circles, xLoc, circleLoc, circleWidth, circleHeight, colorShade):
    i=0
    while i<circles:
        #Ellipse(size 2 array detailing xy, width, height)
        drawEllipse(xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        circleLoc=circleLoc+circleHeight
        i+=1
    return circleLoc
#create a dotgraph, receiving a crosstable and an xLoc
def dotGraph(stuff, xLoc):
    #size is the scale of the ellipse/circle
    size=1
    #circleWidth and circleHeight are the basic sizes of the ellipses (scaled for use in the xy-plane)
    circleWidth=size*.15
    circleHeight=size*0.2
    for c in stuff: 
        #Graph of the first list of data ('case') in our data set
        circleLoc=circleHeight/2
        dictTally=c.genotypeTally()
        try:
            oneCircles= dictTally['0/1']
            ones=True
        except: 
            ones=False
        try:
            twoCircles=dictTally['1/1']
            twos=True
        except:
            twos=False
        if twos:
            colorShade=[0,0,1]
            circleLoc=multiEllipses(twoCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        if ones:
            colorShade=[0,1,0]
            circleLoc=multiEllipses(oneCircles, xLoc, circleLoc, circleWidth, circleHeight, colorShade)
        #plot text above the circles (xcoord, ycoord, text, text placed horizontally left-aligned and 
        #vertically bottom aligned to our location, and having a rotation of 60
        #plt.text(xLoc, circleLoc, c.getpos(), ha='left', va='bottom', rotation=60)
        #Graph of the second list of data ('control') in our data set

        xLoc+=1
    plt.show()
#############################################################

fig=plt.figure()
plt.title("Frequency of Alleles")
plt.xlabel("Chromosome")
ax1 = fig.add_subplot(111)
#plt.axis([0, 10, 0, 10])
ax1.add_line(Line2D([0, 1000], [0, 0],
                  linewidth=1, color='black'))
ax1.grid(True)
fig.subplots_adjust(bottom=0.2)


######################################################
#Begin option parser
parser=OptionParser()
parser.add_option("-f", "--file", dest="filename", help ="data input file", metavar="input")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
#parser.add_option("-c", "--column", dest="column", help="column of data requested", metavar="column")
parser.add_option("-r", "--region", dest="region", help="region requested, chrom: start-stop format", metavar="region")
(options, args)=parser.parse_args()
#Begin regular expression compile for chrom, start, and stop
regionRE=re.compile(r'(.+):(\d+)-(\d+)')

filename=options.filename
filename = './data/458_samples_from_bcm_bi_and_washu.annot.vcf'
options.region="1:157505468-157510468"
m=regionRE.match(options.region)

v = vcf.VCFdata(filename, 1024*1024)
chrom=int(m.groups()[0])
start=int(m.groups()[1])
end=int(m.groups()[2])
chrom=1
start=157505468
end=start+5000
vstuff = v.fetch_range(chrom, start, end)

xlabels=[]
for c in vstuff:
    xlabels.append(c.getpos())
xrange=np.arange(len(xlabels)+1)
plt.xticks(xrange+0.7, xlabels, rotation=65)
y=[-1,0,1]
plt.yticks(y,('case', '','control'), rotation=90)
#insert method to get vector from traits file here
traitstuff=[]
i=1
for v in vstuff:
    crossTable=CrossTable(traitstuff, v.genotypes())
    dotGraph(crossTable, i)
    i=i+0.5

##################################################################################

