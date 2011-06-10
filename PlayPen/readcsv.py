#! /opt/local/bin/python2.6

import rpy2.robjects as robjects
R = robjects.r

from rpy2.robjects.packages import importr 

Rutils = importr('utils')

readcsv = Rutils.read_csv
readtable = Rutils.read_table

fixedfile = "../testing/data/458_traits.fixed"
csvfile = "../testing/data/458_traits.csv"

traitData = readtable(fixedfile, header=True)
traitData2 = readcsv(csvfile, header=True)

for i in range(3):
	print traitData.names[i]
	print traitData[i][0:10]
	print traitData2[i][0:10]

print "=" * 60
# to keep strings as strings rather than factors:
traitData = readtable(fixedfile, header=True)
traitData2 = readcsv(csvfile, header=True)

for i in range(3):
	print traitData.names[i]
	print traitData[i][0:10]
	print traitData2[i][0:10]







