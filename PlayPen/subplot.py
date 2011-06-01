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

filename="domains_LDLR.csv"
column="startcodon"

f=open(filename, "r")
reader=csv.DictReader(f, delimiter=",")
fieldnames=reader.fieldnames

data=[]
listMin=False
listMax=False
for row in reader:
    try:
        newValue = row[column]
        data+=[int(newValue)]
    except:
        print "There is an error with this column.  Your choices are: "
        print fieldnames
        break
    
try:
    data.sort()
    listMin=data[0]
    listMax=data[len(data)-1]
except:
    print "This column could not be graphed.  Perhaps the data are not numbers."
    #set up the bin number and histogram information

#First subplot
#subplot values: first number is the number of rows, second is the number of columns, and third is the number of the graph
subplot(3,1,1)
n, bins, patches = plt.hist(data, len(data)+1, facecolor='green', alpha=0.75, range=(listMin-0.5, listMax+0.5), align='mid', rwidth=1)
plt.xlabel(column)
plt.ylabel('Frequency')
plt.title("Frequency of %s" %column)
plt.grid(True)

#Second subplot
subplot(2,1,2)
n, bins, patches = plt.hist(data, len(data)+1, facecolor='green', alpha=0.75, range=(listMin-0.5, listMax+0.5), align='mid', rwidth=1)
plt.xlabel(column)
plt.ylabel('Frequency')
plt.title("Frequency of %s" %column)
plt.grid(True)

plt.show()
