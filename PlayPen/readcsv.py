#! /opt/local/bin/python2.6

import rpy2.robjects as robjects
R = robjects.r

from rpy2.robjects.packages import importr 

Rutils = importr('utils')

readcsv = Rutils.read_csv
readtable = Rutils.read_table

filename = "../TestData/458_traits.fixed"

traitData = readtable(filename, header=True)


for i in range(3):
	print traitData.names[i]
	print traitData[i][0:10]

print "=" * 60
# to keep strings as strings rather than factors:
traitData = readtable(filename, header=True,stringsAsFactors=False)

for i in range(3):
	print traitData.names[i]
	print traitData[i][0:10]







