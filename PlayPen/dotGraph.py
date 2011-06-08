#!/usr/bin/env python
'''
Created on Jun 8, 2011
A program to generate a dotgraph based on information acquired from a vcf file, using genezoom to create a table of said data
@author: jcc7
'''
import vcf
import CrossTable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

#test data
vecA=('case', 'case', 'control', 'control', 'case', 'control', 'case', 'case', 'case', 'control', 'control')
vecB=('2',      '1',      '0',      '0',      '0',      '1',    '2',    '1',    '2',         '2',        '0')
values=CrossTable.CrossTable(vecA, vecB)
filename = '../TestData/458_samples_from_bcm_bi_and_washu.annot.vcf'
v = vcf.VCFdata(filename, 1024*1024)
pos = 156705468 + 800000
offset = 5000
stuff = v.fetch_range(1, pos, pos + offset)
for c in stuff:
    print (c.getlocus(), c.genotypeTally()) 

##################################################################################
#set up subplots, sharing an x-axis
plt.subplots(nrows=2, ncols=1, sharex=True, frameon=False)
#first subplot setup
plt.subplots_adjust(hspace=0.001)
ax1=plt.subplot(2,1,1)
ax1.axes.get_xaxis().set_visible(False)
ax1.axes.get_yaxis().set_ticks([])
ax1.set_xlim(0,10)
ax1.set_ylim(0,1)
plt.ylabel("Cases")
#second subplot setup
ax2=plt.subplot(2,1,2, sharex=ax1)
ax2.axes.get_yaxis().set_ticks([])
ax2.axes.get_xaxis().set_ticks([])
ax2.set_ylim(-1,0)
plt.ylabel("Controls")
plt.xlabel("Chromosome")

######################################################
#Data graphing
#size is the scale of the ellipse/circle
size=2
#circleWidth and circleHeight are the basic sizes of the ellipses (scaled for use in the xy-plane)
circleWidth=size*.1
circleHeight=size*.03
#Graph of the first list of data ('case') in our data set
circleLoc=circleHeight/2
for a in range(3):
    colorShade=[0,0,a*.3]
    for b in range(values.valueAt(0, a)):
        #Ellipse(size 2 array detailing xy, width, height)
        e=Ellipse([5, circleLoc], circleWidth, circleHeight)
        e.set_facecolor(colorShade)
        ax1.add_artist(e)
        circleLoc=circleLoc+circleHeight
#plot text above the circles (xcoord, ycoord, text, text placed horizontally left-aligned and 
#vertically bottom aligned to our location, and having a rotation of 60
plt.text(5, circleLoc, 'ExampleData',ha='left', va='bottom', rotation=60)
#Graph of the second list of data ('control') in our data set
circleLoc=-circleHeight/2
for a in range(3):
    colorShade=[0,a*.3,0]
    for b in range(values.valueAt(1, a)):
        e=Ellipse([5, circleLoc], circleWidth, circleHeight)
        e.set_facecolor(colorShade)
        ax2.add_artist(e)
        circleLoc=circleLoc-circleHeight

plt.show()