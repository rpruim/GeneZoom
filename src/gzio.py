#!/usr/bin/env python

from gzutils import *
import sys

# Load program constants.
conf_file = find_relative("conf/gz.conf");
execfile(conf_file);


if gz_options['R']:
	try:
		import rpy2.robjects as robjects
		R = robjects.r

		from rpy2.robjects.packages import importr 

		Rutils = importr('utils')
		def post_process_data(data):
			result = {}
			for i in range(len(data)):
				result[data.names[i]] = [ x for x in data[i] ]
			return(result)

		def read_csv(filename, delim=','):
			return post_process_data( Rutils.read_csv(filename, sep=delim, header=True) )

		def read_table(filename, delim=''):
			return post_process_data( Rutils.read_table(filename, sep=delim, header=True) )

	except:
		print sys.stderr >> "Rpy2 unavailable.  Proceeding with alternate tools."
		gzOptions['R'] = False

if not gz_options['R']:
	import numpy as np

	def post_process_data(data):
		result = {}
		for i in range(len(data[0])):
			result[data.dtype.names[i]] = [ x[i] for x in data ]
		return(result)

	def read_csv(filename, delim=','):
		data = np.genfromtxt(filename, delimiter=delim, names=True)
		return post_process_data(data)

	def read_table(filename, delim=' '):
		data = np.genfromtxt(filename, delimiter=delim, names=True)
		return post_process_data(data)

if __name__ == "__main__":
	fixedfile = "../testing/data/458_traits.txt"
	csvfile = "../testing/data/458_traits.csv"

	traitData = read_table(fixedfile)
	traitData2 = read_csv(csvfile)

	for key in traitData:
		print key
		print traitData[key][0:10]

	print "=" * 60
