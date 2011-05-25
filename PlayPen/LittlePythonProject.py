#! /usr/bin/env python 

#author: Jakeniah Christiansen
#A simple program that experiments with OptionParser and making histograms.  May be run as executable, in
#the following form: LittlePythonProject -f nameoffile
#Two text files are available for examples, pivalues and euler, which contain the first few values of pi and euler's number (respectively)
#The program will create a histogram of the density of the digits of these values
from optparse import OptionParser
import csv
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#set up the Option Parser as parser
parser=OptionParser()
parser.add_option("-f", "--file", dest="filename", help ="write report to FILE", metavar="FILE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args)=parser.parse_args()

#read the file into list x

openfile = open(options.filename, 'r')
openfile.seek(0,0)
x=openfile.readlines()


#put list x values into array "values", and store the highest and lowest values
values = []
listMin=0
listMax=0
for i in range(len(x)):
    newValue = x[i]
    values+=[int(newValue)]
    if newValue<listMin:
        listMin=int(newValue)
    if newValue>listMax:
        listMax=int(newValue)

#set up the parameters for the histogram
binNumber=10

# the histogram of the data
n, bins, patches = plt.hist(values, binNumber, normed=0, facecolor='green', alpha=0.75)

plt.xlabel('Digits')
plt.ylabel('Density')
plt.title('Histogram')
xMin=0
xMax=int(binNumber)
yMin=int(listMin)
yMax=int(listMax+1)
plt.axis([xMin, xMax, yMin, yMax])
plt.grid(True)

plt.show()

