#! /usr/bin/env python 

######################################################################################
#
# a program to read in a csv file with headers and produce a histogram of one of the 
# columns.
#
# usage histogram -file somefile.csv -c column_name
#
#author: Jakeniah Christiansen
#date begun: May 25, 2011
######################################################################################

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
parser.add_option("-c", "--column", dest="column", help="write report to FILE")
(options, args)=parser.parse_args()

#initialize the csv reader as reader
#options.filename="domains_LDLR.csv"
#options.column="startcodon"
#reader=csv.reader(open(options.filename))
f=open(options.filename, "r")
dictReader=csv.DictReader(f, delimiter=",")
fieldnames=dictReader.fieldnames

#initialize a list, then input the correspondent value for code for each row in dictReader
data=[]
listMin=False
listMax=False
for row in dictReader:
    newValue = row[options.column]
    data+=[int(newValue)]
    
#sort the list, store the min and max    
data.sort()
listMin=data[0]
listMax=data[len(data)-1]

binNumber=len(data)+1

n, bins, patches = plt.hist(data, binNumber, facecolor='green', alpha=0.75, range=(listMin-0.5, listMax+0.5), align='mid', rwidth=1)

plt.xlabel(options.column)
plt.ylabel('Frequency')
plt.title("Frequency of %s" %options.column)
plt.grid(True)

plt.show()
    
#data=genfromtxt(options.filename, delimiter=",", skip_header=1)
#print data
#for row in data:
 #   for entry in row:
  #      print int(entr

#headers()
#for entry in line:
#    headers.append(entry)
#    print entry   
#i=0
#for entry in headers:
#    print headers[i]
#    i=i+1