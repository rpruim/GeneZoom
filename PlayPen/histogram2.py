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

import csv
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.widgets import RadioButtons

filename="domains_LDLR.csv"
column="startcodon"

f=open(filename, "r")
reader=csv.DictReader(f, delimiter=",")
fieldnames=reader.fieldnames

def dataChange(name):
    
    data=[]
    column=name
    for row in reader:
        try:
            newValue = row[column]
            data+=[int(newValue)]
        except:
            print "There is an error with this column.  Your choices are: "
            print fieldnames
            break
    redraw(data)

def redraw(data):
    #First plot
    ax1.hist(data, facecolor='green', alpha=0.75, align='mid', rwidth=1)
    ax1.grid(True)
    ax1.axes.get_xaxis().set_visible(False)

    #Second plot
    zeros=[0]*len(data)
    ax2.plot(data, zeros, "|")
    ax2.axes.get_yaxis().set_visible(False)


#data=dataChange(column)

#Locations of subplots: left, bottom, width, height (fractions of the whole screen)
rect1=[0.25, 0.15, 0.7, 0.65]
rect2=[0.25, 0.1, 0.7, 0.05]
#sets up the background figure
fig=plt.figure(facecolor='white')
axescolor='#f6f6f6' #(the background color of the axes)
#sets up the location of the two plots as ax1 and ax2
ax1=fig.add_axes(rect1, axisbg=axescolor)
ax2=fig.add_axes(rect2, axisbg=axescolor, sharex=ax1)

#button control
controlPanel=fig.add_axes([0.05, 0.15, 0.15, 0.25])
radio=RadioButtons(controlPanel, fieldnames, active=0)
radio.on_clicked(dataChange)
dataChange(column)

show()
