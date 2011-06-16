#!/ usr/bin/env python
#
# utilities for interacting with the UCSC Genome Browser data.
# based on a similar set of utilities written in perl by Peter Chines.
#

#use Carp qw(croak confess);

import re
import sys
from gzio import *
from gzutils import *
import logging

class BED:
	def __init__(self, filename, keys=None):
		fh = multiopen( filename )
		self._keys = keys
		self._key2index = {}
		index = 0
		try:
			for key in keys:
				self._key2index[key] = index
				index = index + 1
		except TypeError:
			self._key2index = {}

		self._rows = [BEDrow(s, parent=self) for s in fh.readlines() ]

	def __getitem__(self, r, c=None):
		if c == None:
			return self._rows[r]
		return self._rows[r][c]

	def get_all_rows(self):
		return [r for r in self._rows ]

	def get_rows(self, name='', exact=True):
		if name == '':
			exact=False
		if exact:
			return [r for r in self._rows if name == r['name'] ]
		else:
			return [r for r in self._rows if re.search(name, r['name']) ]

class BEDrow:
	def __init__(self, s, parent=None):
		self._contents = s.strip().strip(',').split('\t') 
		self._parent = parent

	def __repr__(self):
		return ":".join(self._contents)

	def __getitem__(self, item):
		if type(item) in [int, long]:
			return self._contents[item]
		return self._contents[self._parent._key2index[item]]
	

	def get_txRange(self):
		return (self['txStart'], self['txEnd'])

	def get_utr(self):
		result = [ (min(a,b), max(a,b)) for a, b in zip(self.get_txRange(), self.get_cdsRange()) ]
		for i in range(2):
			if result[i][0] == result[i][1]:
				result[i] = tuple()
		return result

	def get_cdsRange(self):
		return (self['cdsStart'], self['cdsEnd'])

	def get_exonStarts(self):
		return [int(a) for a in self['exonStarts'].split(',') if len(a) > 0]

	def get_exonEnds(self):
		return [int(a) for a in self['exonEnds'].split(',') if len(a) > 0]

	def get_exons(self):
		return [ (a,b) for a,b in zip( self.get_exonStarts(), self.get_exonEnds() ) ]


if __name__ == "__main__":

	refFlatKeys = ['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
	knownGeneKeys = ['name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds',
			'proteinID','alignID']
	refFlat = BED('../testing/data/refFlat.txt.gz.1', keys=refFlatKeys)
	knownGene = BED('../testing/data/knownGene.txt', keys=knownGeneKeys)

	for i in range(3):
		r = refFlat[i]
		print r['name'], r['chrom'], r['strand'] 
		print '\t txRange:  ',r.get_txRange() 
		print '\tcdsRange:  ',r.get_cdsRange()
		print '\t   exons: ', r.get_exons() 
		print '\t     UTR: ', r.get_utr()
		print '\n'
	print refFlat.get_rows('NR_026820')
	print "=" * 60
	for i in range(3):
		r = knownGene[i]
		print r['name'], r['chrom'], r['strand'] 
		print '\t txRange:  ',r.get_txRange() 
		print '\tcdsRange:  ',r.get_cdsRange()
		print '\t   exons: ', r.get_exons() 
		print '\t     UTR: ', r.get_utr()
		print '\n'

	print knownGene.get_rows('NR_026820')
	
