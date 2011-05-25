
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
(options, args)=parser.parse_args()

#initialize the csv reader as reader
options.filename="domains_LDLR.csv"
#reader=csv.reader(open(options.filename))
f=open(options.filename, "r")
dictReader=csv.DictReader(f, delimiter=",")
fieldnames=dictReader.fieldnames

#initialize a list, then input the correspondent value for code for each row in dictReader
code=[]
for row in dictReader:
    newValue = row['code']
    code+=[int(newValue)]
#sort the list, store the min and max    
#code.sort()
#listMin=code[0]
#listMax=code[len(code)-1]

binNumber=3

n, bins, patches = plt.hist(code, binNumber, normed=0, facecolor='green', alpha=0.75, align='mid', rwidth=1, range=(1, 3))

plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title("Code Frequency")
xMin=0.5
xMax=3.5
yMin=0
yMax=15
plt.axis([xMin, xMax, yMin, yMax])
plt.grid(True)

plt.show()
#print fieldnames
#for row in dictReader:
#    print row['endcd']
    
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