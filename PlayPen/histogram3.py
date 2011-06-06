#! /usr/bin/env python 

######################################################################################
#
# A program for practice combining histograms, subplots, and tick marks
# (practice using matplotlib in general)
#
# usage histogram -file somefile.csv -c column_name
#
#author: Jakeniah Christiansen
#date of creation: June 1, 2011
######################################################################################
from optparse import OptionParser
import csv
import matplotlib.pylab as pylab
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import Slider, show, draw, axes
from matplotlib.widgets import RadioButtons

class Histogram:
    #Initial setup. initialize reader, fieldnames, and prepare data
    def __init__(self, filename, column, format):
        #if the file is a csv, then initialize the csv reader as reader
        if format=="csv":
            f=open(filename, "r")
            self.reader=csv.DictReader(f, delimiter=",")
            self.fieldnames=self.reader.fieldnames
        elif format=="vcf":
            #Here is where the code for initiliazing a vcf reader would go
            print "This works (kinda)."
        else:
            print "Not a valid file type.  Please use either csv or vcf."
        self.column=column
        self.dataSetup(column)
    
    def dataSetup(self, column):  
        self.data=[]
        for row in self.reader:
            newValue = row[column]
            try:
                if is_number(newValue):
                    self.data+=[int(newValue)]
                else:
                    self.data+=[newValue]
            except:
                print "There is an error with this column.  Your choices are: "
                print self.fieldnames
                break
        self.data.sort()
        self.min=self.data[0]
        self.max=self.data[len(self.data)-1]

    def length(self):
        return len(self.data)
    def min(self):
        return self.min
    def max(self):
        return self.max
    def data(self):
        return self.data
    def fieldnames(self):
        return self.fieldnames
    def column(self):
        return self.column
    
######################################################################################    
#A method to check whether or not a value is a number (to see if it is histogram graphable or not)
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#set up the Option Parser as parser
parser=OptionParser()
parser.add_option("-f", "--file", dest="filename", help ="write report to FILE", metavar="FILE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option("-c", "--column", dest="column", help="write report to FILE")
#parser.add_option("-l", "--location", dest="location", help="center of requested data")
(options, args)=parser.parse_args()

#Take these lines out later:
filename="domains_LDLR.csv"
column="startcodon"
format="csv"
location="500"
    
#if options.filename.endswith(".csv"):
#    format="csv"
#elif options.filename.endswith(".vcf"):
#    format = "vcf"
#else:
#    format="invalid"



#Histogram constructor
hist=Histogram(filename, column, format)

#Setup the position of the plots (left, bottom, width, height)
histRect=[0.15, 0.25, 0.7, 0.65]
rugRect=[0.15, 0.20, 0.7, 0.05]
#controlRect=[0.05, 0.20, 0.15, 0.25]
locRect=[0.15, 0.1, 0.65, 0.03]
zoomRect=[0.15, 0.05, 0.65, 0.03]

fig=plt.figure(facecolor='white')
axescolor="#f6f6f6"
#controlPanel=fig.add_axes(controlRect)
#controlFields=[]
#radio=RadioButtons(controlPanel, hist.fieldnames, active=0)
axloc  = axes(locRect, axisbg=axescolor)
axzoom = axes(zoomRect, axisbg=axescolor)
initLoc=int(location)
#initLoc=(hist.max-hist.min)/2+hist.min
initWidth=500
sLoc = Slider(axloc, 'Location', hist.min, hist.max , valinit=initLoc)
sZoom = Slider(axzoom, 'Width', 1, 1000, valinit=initWidth)

#set up the axes
ax1=fig.add_axes(histRect, axisbg=axescolor)
ax1.grid(True)
plt.title("Frequency of %s" %column)
plt.ylabel('Frequency')
ax1.axes.get_xaxis().set_visible(False)
ax2=fig.add_axes(rugRect, axisbg=axescolor, sharex=ax1)
zeros=[0]*hist.length()
ax2.axes.get_yaxis().set_visible(False)
plt.xlabel(column)

#draw the initial histograms
ax1.hist(hist.data, facecolor='green', alpha=0.75, align='mid', range=(hist.min, hist.max), rwidth=1)
ax2.plot(hist.data, zeros, "|")

#if an alteration in the requested slider values occurs, then create the new histogram based on these values
def update(val):
    ax1.clear()
    ax2.clear()
    xLoc = sLoc.val
    xWidth = sZoom.val
    leftValue=xLoc-(xWidth/2)
    rightValue=xLoc+(xWidth/2)
    ax1.hist(hist.data, facecolor='green', alpha=0.75, align='mid', range=(leftValue, rightValue), rwidth=1)
    ax2.plot(hist.data, zeros, "|")
    draw()
sLoc.on_changed(update)
sZoom.on_changed(update)

show()
#save to a pdf file
fig.savefig('Example Histogram-Frequency of %s.pdf' % column)
