'''
Created on July 11, 2011
A module for reading compressed, tab-delimited files using tabix-created tbi files.
@author: rpruim
'''

import string
import zlib
from gzio import *
import re
import struct
import logging
import vcfrow

def readChar(fh):
	return struct.unpack( '<c', fh.read(1) )[0]

def readInt8(fh):
	return struct.unpack( '<b', fh.read(1) )[0]

def readUInt8(fh):
	return struct.unpack( '<B', fh.read(1) )[0]

def readInt16(fh):
	return struct.unpack( '<h', fh.read(2) )[0]

def readUInt16(fh):
	return struct.unpack( '<H', fh.read(2) )[0]

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


class tabixReader:

	def __init__(self, filename):
		self._bgzf_filename = filename
		self._tbi_filename = filename + '.tbi'
		self.open()

	def open(self):
		try:
			self._bgzf = open(self._bgzf_filename)
			self._tbi = multiopen(self._tbi_filename)
			self._magic = self._tbi.read(4)
			assert( self._magic == 'TBI\x01' )
			self._n_ref = readInt32(self._tbi)
			self._format = readInt32(self._tbi)
			self._col_seq = readInt32(self._tbi)
			self._col_beg = readInt32(self._tbi)
			self._col_end = readInt32(self._tbi)
			self._meta = readInt32(self._tbi) # chr( struct.unpack('<i', self._tbi.read(4) )[0] )
			self._skip = readInt32(self._tbi)
			self._l_nm = readInt32(self._tbi)
			format = '<' + str(self._l_nm) + 's'
			self._seq_names = struct.unpack(format, self._tbi.read(self._l_nm) )[0]
			# names are 0 terminated, so split on \x00 and drop last empty name
			self._seq_names = self._seq_names.split('\x00')[:-1]    
		except Exception as e:
			msg = 'invalid tbi file [' + repr(e) + ']'
			die(msg)
		self._index= {}

		for name in self._seq_names:
			n_bin = readInt32(self._tbi)
			self._index[name] = {'n_bin': n_bin, 'bins':[], 'bin': {}, 'intervals': [] }
			for b in range(n_bin):
				bin_nr = readUInt32(self._tbi)
				n_chunk = readInt32(self._tbi)
				for c in range(n_chunk):
					cnk_beg = readUInt64(self._tbi)
					cnk_end = readUInt64(self._tbi)
					self._index[name]['bins'].append( 
						(bin_nr, c, n_chunk, cnk_beg >> 16, cnk_beg % (1<<16), cnk_end >> 16, cnk_end % (1<<16), cnk_end - cnk_beg) 
					)
					self._index[name]['bin'][bin_nr] = \
						{ 'nr': bin_nr,
						  'chunk': c,
						  'num_chunks': n_chunk,
						  'chunk_begin_coffset': cnk_beg >> 16,         # high order bits point to beginning of bin
						  'chunk_begin_uoffset': cnk_beg % (1<<16),     # lo order bits locate chunk within bin
						  'chunk_end_coffset': cnk_end >> 16,         # high order bits point to beginning of bin where chunk ends
						  'chunk_end_uoffset': cnk_end % (1<<16)      # lo order bits locate chunk end with bin
						} 
			n_intv = readInt32(self._tbi)
			for i in range(n_intv):
				ioff = readUInt64(self._tbi)
				self._index[name]['intervals'].append((i,ioff)) 

		self._tbi.close()

	def close(self):
		self._tbi.close()
		self._bgzf.close()

	def get_tbi_header(self):

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
				'names': self._seq_names
				}

	def bin2offset(self, chrom, bin):
		try:
			return	self._index[chrom]['bin'][bin]
		except KeyError:
			return None

	def coffset2string(self, coffset):
		self._bgzf.seek( coffset ) # bin_info['chunk_begin_coffset'] )
		try:
			assert( readUInt8(self._bgzf) == 31 )
			assert( readUInt8(self._bgzf) == 139 )
			assert( readUInt8(self._bgzf) == 8 )
			assert( readUInt8(self._bgzf) == 4 )
			mod_time = readUInt32(self._bgzf)
			xfl = readUInt8(self._bgzf)
			os = readUInt8(self._bgzf)
			xlen = readUInt16(self._bgzf)
			assert( readUInt8(self._bgzf) == 66 )
			assert( readUInt8(self._bgzf) == 67 )
			assert( readUInt16(self._bgzf) == 2 )
			bsize = readUInt16(self._bgzf)  # total block size minus 1

		except AssertionError as e:
			logging.critical( "Cannot resolve bgzf file format issues.")
			logging.critical( repr(e) )
		self._bgzf.seek( coffset ) # bin_info['chunk_begin_coffset'] )
		s = zlib.decompress( self._bgzf.read(  bsize + 1), zlib.MAX_WBITS + 32 )
		return s

	def bin2string(self, chrom, bin):
		bin_info = self.bin2offset(chrom, bin)
		if bin_info == None:
			return None

		s = self.coffset2string(bin_info['chunk_begin_coffset'])
		if  bin_info['chunk_begin_coffset'] == bin_info['chunk_end_coffset']:
			return s[ bin_info['chunk_begin_uoffset']:bin_info['chunk_end_uoffset'] ]

		logging.debug('Should add a consistency check here when spanning multiple blocks.')
		t = self.coffset2string(bin_info['chunk_end_coffset'])
		return s[ bin_info['chunk_begin_uoffset']: ] + t[ :bin_info['chunk_end_uoffset'] ] 

	def bin2vcf(self, chrom, bin):
		s = self.bin2string(chrom, bin)
		if s == None:
			return []
		return [ vcfrow.VCFrow(x) for x in s.split('\n') ] 

	'''
	returns a list of strings, one per bin
	'''
	def reg2strings(self, chrom, start, end):
		bins = reg2overlappingBins(start, end)
		result = []
		for b in bins:
			s = self.bin2string(chrom, b)
			if s != None:
				result.extend( s.split('\n') )
		return result
		
	def reg2list(self, chrom, start, end):
		return [ s.split('\t') for s in self.reg2strings(chrom, start, end) ]

	def reg2vcf(self, chrom, start, end):
		return[ vcfrow.VCFrow(s) for s in self.reg2strings(chrom,start,end) ]

