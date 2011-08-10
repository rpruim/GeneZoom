#!/usr/bin/env python

from gzutils import *
import sys

# Load program constants.
conf_file = find_relative("conf/gz.conf");
if conf_file:
	execfile(conf_file)
else:
	warn('Unable to find configuration file.  Proceeding anyway.')
	gz_options = {}

gz_options['R'] = False

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
			return post_process_data( Rutils.read_csv(filename, sep=delim, header=True, stringsAsFactors=False) )

		def read_table(filename, delim=''):
			return post_process_data( Rutils.read_table(filename, sep=delim, header=True, stringsAsFactors=False) )

	except:
		print sys.stderr >> "Rpy2 unavailable.  Proceeding with alternate tools."
		gzOptions['R'] = False

if not gz_options['R']:
	try:
		#assert (False == True)
		import numpy as np
		print "Using numpy"

		def post_process_data(data):
			result = {}
			for i in range(len(data[0])):
				result[data.dtype.names[i]] = [ x[i] for x in data ]
			return(result)

		def post_process_data(data):
			result = {}
			for i in range(len(data[0])):
				result[data.dtype.names[i]] = [ x[i] for x in data ]
			return(result)

		def read_csv(filename, delim=','):
			data = np.genfromtxt(filename, delimiter=delim, names=True, dtype=None)
			return post_process_data(data)

		def read_table(filename, delim=' '):
			data = np.genfromtxt(filename, delimiter=delim, names=True, dtype=None)
			return post_process_data(data)

	except:
		import csv as csv

		def post_process_data(data, headers):
			result = {}
			for i in range(len(headers)):
				try: 
					result[h] =  [ int(x[i]) for x in data ]
				except:
					try:
						result[headers[i]] =  [ float(x[i]) for x in data ]
					except:
						result[headers[i]] =  [ x[i] for x in data ]

			return(result)

		def read_csv(filename, delim=','):
			fh = open(filename,'r')
			dialect = csv.Sniffer().sniff(fh.read(1024)) 
			fh.seek(0)
			if csv.Sniffer().has_header(fh.read(1024)):
				fh.seek(0)
				reader = csv.reader(fh, dialect) 
				headers = reader.next()
			else:
				logging.warning('Unable to interpret csv header row. Proceeding with generic names.')
				fh.seek(0)
				reader = csv.reader(fh, dialect) 
				headers = [ "V"+str(i) for i in range(len(data[0])) ]
			data = [row for row in reader]
			return post_process_data(data, headers)

		def read_table(filename, delim=' '):
			return read_csv(filename, delim=delim)

def multiopen( filename, format='r' ):
	try:
		import gzip
		filehandle = gzip.open(filename, format)                  
		filehandle.read(1024)
		filehandle.seek(0)
	except IOError, ImportError:
		try:
			import bz2
			filehandle = bz2.BZ2File(filename, format)                  
			filehandle.read(1024)
			filehandle.seek(0)
		except IOError, ImportError:
			filehandle = open(filename, format)                  
			filehandle.read(1024)
			filehandle.seek(0)

	return(filehandle)

if __name__ == "__main__":
	fixedfile = "../testing/data/458_traits.txt"
	csvfile = "../testing/data/458_traits.csv"

	traitData = read_table(fixedfile)
	traitData2 = read_csv(csvfile)

	for key in traitData:
		print key
		print traitData[key][0:10]
		assert( traitData[key][0:10] == traitData2[key][0:10] )

	print "=" * 60
