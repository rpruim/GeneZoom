###
### tools for dealing with vcf data
###

import tabix
import string
from gzio import *
import re

class vcfRow:
	def __init__(self, list):
#		s.strip()
		self._items = list  # s.split('\t')

	def get_chrom(self):
		return self._items[0]

	def get_pos(self):
		return int(self._items[1])

	def get_locus(self):
		return ( int(self._items[0]),  int(self._items[1]) )

	def get_genotypes(self, start=9, stop=None):
		if stop == None:
			stop = len(self._items)
		return [ string.split(g,":")[0] for g in self._items[start:stop] ]

	def genotypeTally(self, start=9, stop=None):
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
		return repr(self._items)   #string.join(self._items, "\t")


class vcfReader:
	headerRE     = re.compile(r"(.+?)(#CHROM[^\n]+)", re.M + re.DOTALL)
	
	def __init__(self, filename):
		
		self._filename = filename
		self._blocksize = 1024 * 1024

		# grab header and meta information
		attempts = 0
		self._filehandle = multiopen(self._filename)
		while attempts < 10:
			attempts = attempts + 1
			self._filehandle.seek(0)
			top_of_file = self._filehandle.read(self._blocksize)
			m = vcfReader.headerRE.search(top_of_file)
			if not m:
				self._blocksize = self._blocksize * 4
				continue
			self._metaStr = m.group(1)
			self._headerStr = m.group(2)
			self._headers = self._headerStr.split('\t')
			self._meta = self._metaStr.split('\n')
			self._num_cols = len(self._headers)
			break
		self._filehandle.close()
		self._vcf = tabix.Tabix(filename)

	def close(self):
		self._filehandle.close()

	def get_headers(self):
		return self._headers

	def get_metaStr(self):
		return(self._metaStr)

	def get_ids(self):
		return self._headers[9:]

	def query(self, chrom, start, stop):
		return [ vcfRow(row) for row in self._vcf.query(chrom, start, stop) ]

	def raw_query(self, chrom, start, stop):
		return self._vcf.query(chrom, start, stop) 

def extract_genotypes(row):
	return [ a.split(":")[0] for a in row[9:] ]

def extract_chrom(row):
	return row[0]

def extract_pos(row):
	return row[1]

def tally(x):
	result = {}
	for item in x:
		try:
			result[ item ] = 1 + result[ item ]
		except KeyError:
			result[ item ] = 1 

	return result

def xtally(x,y):
	result = {}
	for item in zip(x,y):
		try:
			result[ item ] = 1 + result[ item ]
		except KeyError:
			result[ item ] = 1 

	return result


if __name__ == '__main__':
	import tabix
	import string

	filename = '../testing/data/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf.gz'
	filename = '../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1'


	v = vcfReader(filename)
#	print v.get_ids()
	print len(v.get_ids()), " subjects"

#	v = tabix.Tabix(filename)

	if True:
		pos = 156705468 + 800000
		offset = 5000

		if False:
			print 'query(', pos, pos+offset, ')...'
			markers = [locus for locus in v.raw_query('1', pos, pos + offset)]
			for c in markers:
				print  extract_chrom(c), extract_pos(c), ": ",
				print  tally( extract_genotypes(c) )
			print '\nTotal: ', len(markers), 'markers'

		else:
			print 'query(', pos, pos+offset, ')...'
			markers = v.query('1', pos, pos + offset)
			for c in markers:
				print c.get_locus(),
				print tally( c.get_genotypes() )
				print c.get_locus(),
				print c.genotypeTally() 
			print '\nTotal: ', len(markers), 'markers'


	filename = '../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1'
	v = vcfReader(filename)
	check( v.get_ids()[0] == 'HG00096', "Checking vcfReader.get_ids()") 
	check( len(v.get_ids()) == 458 )
	pos = 156705468 + 800000
	offset = 5000
	markers = v.query('1', pos, pos + offset)
	check( len(markers) == 20 , "Checking vcfReader.query()" )
	check( markers[0].genotypeTally()['0/0'] == 155 , "Checking vcfRow.genotypeTally()")
	check( markers[0].get_locus() == (1, 157508826), "Checking vcfRow.get_locus()") 
	check( markers[0].get_pos() == 157508826, "Checking vcfRow.get_pos()") 
	check( markers[0].get_chrom() == '1', "Checking vcfRow.get_chrom()") 
	print "\nDone."
