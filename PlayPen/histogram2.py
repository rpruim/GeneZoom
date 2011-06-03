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
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import *
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

(options, args)=parser.parse_args()

#Take these lines out later:
filename="domains_LDLR.csv"
column="startcodon"
format="csv"
#if options.filename.endswith(".csv"):
#    format="csv"
#elif options.filename.endswith(".vcf"):
#    format = "vcf"
#else:
#    format="invalid"

#Histogram constructor
hist1=Histogram(filename, column, format)

#Setup for the graph (left, bottom, width, height0
histRect=[0.25, 0.15, 0.7, 0.65]
rugRect=[0.25, 0.1, 0.7, 0.05]
controlRect=[0.05, 0.15, 0.15, 0.25]
fig=plt.figure(facecolor='white')
axescolor="#f6f6f6"
controlPanel=fig.add_axes(controlRect)
controlFields=[]
radio=RadioButtons(controlPanel, hist1.fieldnames, active=0)


ax1=fig.add_axes(histRect, axisbg=axescolor)
ax1.grid(True)
ax1.axes.get_xaxis().set_visible(False)
ax2=fig.add_axes(rugRect, axisbg=axescolor, sharex=ax1)
zeros=[0]*hist1.length()
ax2.axes.get_yaxis().set_visible(False)
ax1.hist(hist1.data, facecolor='green', alpha=0.75, align='mid', range=(hist1.min, hist1.max), rwidth=1)
ax2.plot(hist1.data, zeros, "|")

def columnChange(column):

    hist1=Histogram(filename, column, format)
    ax1.clear()
    ax2.clear()
    ax1.hist(hist1.data, facecolor='green', alpha=0.75, align='mid', range=(hist1.min, hist1.max), rwidth=1)
    ax2.plot(hist1.data, zeros, "|")
    
radio.on_clicked(columnChange)
show()
