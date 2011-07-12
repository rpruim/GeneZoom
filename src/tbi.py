'''
Created on July 11, 2011
A module for reading tabix-created tbi files.
@author: rpruim
'''

import string
from gzio import *
import re
import struct

def readInt32(fh):
	return struct.unpack( '<i', fh.read(4) )[0]

def readUInt32(fh):
	return struct.unpack( '<I', fh.read(4) )[0]

def readInt64(fh):
	return struct.unpack( '<q', fh.read(8) )[0]

def readUInt64(fh):
	return struct.unpack( '<Q', fh.read(8) )[0]

def reg2bin(beg, end):
	end = end - 1
	if (beg >> 14) == (end >> 14): return (( 1<<15)-1)/7 + (beg>>14)
	if (beg >> 17) == (end >> 17): return (( 1<<12)-1)/7 + (beg>>17)
	if (beg >> 20) == (end >> 20): return (( 1<<9)-1)/7 + (beg>>20)
	if (beg >> 23) == (end >> 23): return (( 1<<6)-1)/7 + (beg>>23)
	if (beg >> 26) == (end >> 26): return (( 1<<3)-1)/7 + (beg>>26)
	return(0)

def reg2overlappingBins(rbeg, rend):
	rend = rend - 1
	result = [0]
	for k in range( 1 + (rbeg>>26), 1 + 1 + (rend>>26) ):
		result.append(k)
	for k in range( 9 + (rbeg>>23), 1 + 9 + (rend>>23) ):
		result.append(k)
	for k in range( 73 + (rbeg>>20), 1 + 73 + (rend>>20) ):
		result.append(k)
	for k in range( 585 + (rbeg>>17), 1 + 585 + (rend>>17) ):
		result.append(k)
	for k in range( 4681 + (rbeg>>14), 1 + 4681 + (rend>>14) ):
		result.append(k)
	return(result)
class tbiReader:


	def __init__(self, filename):
		
		self._filename = filename
		self._filehandle = multiopen(self._filename)
		self._magic = self._filehandle.read(4)
		self._n_ref = readInt32(self._filehandle)
		self._format = readInt32(self._filehandle)
		self._col_seq = readInt32(self._filehandle)
		self._col_beg = readInt32(self._filehandle)
		self._col_end = readInt32(self._filehandle)
		self._meta = chr( struct.unpack('<i', self._filehandle.read(4) )[0] )
		self._skip = readInt32(self._filehandle)
		self._l_nm = readInt32(self._filehandle)
		format = '<' + str(self._l_nm) + 's'
		self._names = struct.unpack(format, self._filehandle.read(self._l_nm) )[0]
		# names are 0 terminated, so split on \x00 and drop last empty name
		self._names = self._names.split('\x00')[:-1]    
		self._data= {}

		for name in self._names:
			n_bin = readInt32(self._filehandle)
			self._data[name] = {'n_bin': n_bin, 'bins':[], 'bin': {}, 'intervals': [] }
			for b in range(n_bin):
				bin_nr = readUInt32(self._filehandle)
				n_chunk = readInt32(self._filehandle)
				for c in range(n_chunk):
					cnk_beg = readUInt64(self._filehandle)
					cnk_end = readUInt64(self._filehandle)
					self._data[name]['bins'].append( (bin_nr, n_chunk, cnk_beg, cnk_end) )
					self._data[name]['bin'][bin_nr] =  (n_chunk, cnk_beg, cnk_end) 
			n_intv = readInt32(self._filehandle)
			for i in range(n_intv):
				ioff = readUInt64(self._filehandle)
				self._data[name]['intervals'].append(ioff) 

	def close(self):
		self._filehandle.close()

	def get_header_info(self):

		return { 
				'magic': self._magic,
				'n_ref': self._n_ref,
				'col_seq': self._col_seq,
				'col_beg': self._col_beg,
				'col_end': self._col_end,
				'meta': self._meta,
				'skip': self._skip,
				'format': self._format,
				'l_nm': self._l_nm,
				'names': self._names
				}


	def reg2offset(self, chrom, start, end):
		bins = reg2bins(start, end)
		return [ (a,b,c,d) for (a,b,c,d) in self._data[chrom]['bins'] if a in reg2bins(start,end) ]


