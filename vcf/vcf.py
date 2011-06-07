###
### Version 0.1b
### vcf.py - a prototype Python API for VCF parsing
###
### Written by Mikhail Spivakov 2011, current version distributed under the terms of CC BY-SA licence (Creative Commons Attribution-ShareAlike)
### Please report bugs/comments/suggestions to spivakov@ebi.ac.uk
### 
###

import re
import os
import sys
import cPickle
import random
import pprint
import linecache
import string
import gzip
import bz2

# Creates subdictionaries from dictionaries containing only a specified subset of the original keys  
extract = lambda keys, dict: reduce(lambda x, y: x.update({y[0]:y[1]}) or x, map(None, keys, map(dict.get, keys)), {})

# counts number of lines in a file
def bufcount(filename, max=100000):
    f = open(filename)                  
    lines = 0
    buf_size = 64 * 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

def multiopen( filename, format ):
	try:
		filehandle = gzip.open(filename, format)                  
		filehandle.read(1024)
		filehandle.seek(0)
	except IOError:
		try:
			filehandle = bz2.BZ2File(filename, format)                  
			filehandle.read(1024)
			filehandle.seek(0)
		except IOError:
			filehandle = open(filename, format)                  
			filehandle.read(1024)
			filehandle.seek(0)

	return(filehandle)

class VCFdata:
	headerRE     = re.compile(r'#CHROM')
	headerRE     = re.compile(r"(?P<meta>.+)(?P<header>#CHROM[^\n]+)\n", re.M)
	headerRE     = re.compile(r"(.+?)(#CHROM[^\n]+)", re.M + re.DOTALL)
	firstLocusRE = re.compile(r'\n(\d+)\t([^\t]+)\t')

	def __init__ (self, filename, blocksize=1024*1024):
		self._filename = filename
		self._filesize = os.path.getsize(filename)
		self._filehandle = multiopen(filename,'r')

		self._blocksize = blocksize # should be large enough to contain meta and header info; stretched if not.
		self._headerStr = None
		self._metaStr = None
		self._read_from_file = self._filehandle.read

		# grab header and meta information
		attempts = 0
		while attempts < 10:
			attempts = attempts + 1
			self._filehandle.seek(0)
			top_of_file = self._read_from_file(self._blocksize)
			m = VCFdata.headerRE.search(top_of_file)
			if not m:
				self._blocksize = self._blocksize * 4
				continue
			self._metaStr = m.group(1)
			self._headerStr = m.group(2)
			self._headers = self._headerStr.split('\t')
			self._num_cols = len(self._headers)
			self._data_start = len(self._metaStr) + len(self._headerStr)
			self._num_blocks = self._filesize / self._blocksize
			break

	def fetch_block(self, block, raw=False):
		self._filehandle.seek(self._data_start + block * self._blocksize)
		blockStr = self._read_from_file(self._blocksize)
		if raw:
			return blockStr
		if len(blockStr) < 1:
			return []
		while not blockStr[-1] == '\n':
			blockStr = blockStr + self._read_from_file(1)
		# block = [ l.split('\t') for l in blockStr.split('\n') ]
		block = [ VCFrow(l) for l in blockStr.split('\n') ]
		block = [ r for r in block if len(r) == self._num_cols ]
		return block

	# optimize this later by just grabbing position?
	def first_locus_in_block(self, block):
		try:
			blockStr = self.fetch_block(block, raw=True)
			m = VCFdata.firstLocusRE.search(blockStr)
			chrom = int( m.groups()[0] )
			pos = int( m.groups()[1] )
			return ( chrom, pos )
		except Exception as e: 
			print >> sys.stderr, "Trouble in first_locus_in_block(", block, ")"
			print >> sys.stderr, e
			print >> sys.stderr, type(e)
			return None

	def first_locus_in_block_alt(self, block):
		self._filehandle.seek(self._data_start + block * self._blocksize)
		while not self._filehandle.read(1) == '\n':
			pass
		chunk = self._filehandle.read(100)
		firstLocusRE = re.compile(r'([^\t]+)\t([^\t]+)\t')
		m = firstLocusRE.search(chunk)
		chrom = int( m.groups()[0] )
		pos = int( m.groups()[1] )
		return ( chrom, pos )

	def last_locus_in_block(self, block):
		try:
			block = self.fetch_block(block)
			return block[-1].getlocus()
		except IndexError as e:
			print >> sys.stderr, "Trouble in last_locus_in_block(", block, ")"
			print >> sys.stderr, e
			print >> sys.stderr, type(e)
			return None

	def fetch_range(self, chrom, start, stop):
		search_space = [0, self._num_blocks-1]
		start_locus = (chrom, start)
		stop_locus = (chrom, stop)
		while search_space[0] + 1 < search_space[1]:
			#print >> sys.stderr, " ", str(search_space), "." ,
			mid = int( sum(search_space) / 2 )
			#print >> sys.stderr,  ".." ,
			if self.first_locus_in_block(mid) > start_locus:
				search_space[1] = mid
			else:
				search_space[0] = mid
			#print >> sys.stderr, "...",
		#print >> sys.stderr, "  "
		start_block = search_space[0]

		search_space = [start_block, self._num_blocks-1]
		while search_space[0] + 1 < search_space[1]:
			#print >> sys.stderr, "  ", str(search_space),
			mid = int( sum(search_space) / 2 )
			if self.first_locus_in_block(mid+1) < stop_locus:
				search_space[0] = mid
			else:
				search_space[1] = mid
		#print >> sys.stderr, "  "

		stop_block = search_space[0]

		results = []
		for block in range(start_block, 1+stop_block):
			results.extend( self.fetch_block(block) )

		return [r for r in results if r.getlocus() >= start_locus and r.getlocus() <= stop_locus]

	# return a string summarizing the data contained in VCF
	def __repr__(self):
		repr = "VCF file"
#		repr = repr + self._headerStr
		repr = repr + '\n\t' + str(self._filesize / (1024 * 1024) ) + ' Mb'
		repr = repr + '\n\t' + str(self._num_blocks) + ' blocks'
		return repr

	def select(self, chrom, start, stop):
		return 'Size: ' + str(self._filesize)

class VCFrow:
	def __init__(self, s):
		s.strip()
		self._items = s.split('\t')

	def getchrom(self):
		return self._items[0]

	def getpos(self):
		return self._items[1]

	def getlocus(self):
		return ( int(self._items[0]),  int(self._items[1]) )

	def genotypes(self, start=9, stop=None):
		if stop == None:
			stop = len(self._items)
		result = {}
		for col in range(start, stop):
			try:
				keys = string.split(self._items[col], ":")
				result[ keys[0] ] = 1 + result[ keys[0] ]
			except KeyError:
				result[ keys[0] ] = 1 

		return result

	def __len__(self):
		return len(self._items)

	def __repr__(self):
		return string.join(self._items, "\t")

def coord (c):
        ccoord = re.search("(.+)@(.+)", c)
        return (ccoord.group(1), ccoord.group(2))

if __name__ == '__main__':
	filename = '../TestData/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf'
	filename = '../TestData/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf.gz'
	filename = '../TestData/458_samples_from_bcm_bi_and_washu.annot.vcf.gz'
	filename = '../TestData/458_samples_from_bcm_bi_and_washu.annot.vcf'
#	v = VCFdata(filename, 1024*64)
	v = VCFdata(filename, 1024*1024)
	print v
#	print v.first_locus_in_block(0)
#	print v.last_locus_in_block(0)

#	print v.first_locus_in_block(1)
#	print v.last_locus_in_block(1)

#	print v.first_locus_in_block(5649)
#	print v.last_locus_in_block(5649)

#	bl = 4
#	print 'block ' + str(bl), 'contains', 
#	contents = v.fetch_block(bl)
#	print len( contents ), # [ c.getlocus() for c in contents ] 
#	print "loci."

#	print "\tFirst Locus:  ",  v.first_locus_in_block(bl)

	if False:
		for i in range(v._num_blocks):
			print >> sys.stderr, i, v.first_locus_in_block(i)

	if True:
		pos = 156705468 + 800000
		offset = 5000
		stuff = v.fetch_range(1, pos, pos + offset)
		print 'fetch_range(', pos, pos+offset, ')...'
		# print string.join( [str(l) for l in stuff], "\n" )
		for c in stuff:
			print (c.getlocus(), c.genotypes()) 
		print '\nTotal: ', len(stuff), 'markers'


