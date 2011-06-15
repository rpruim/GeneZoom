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
	def __init__(self, filename):
		fh = multiopen( filename )
		self._rows = [BEDrow(s) for s in fh.readlines() ]

	def __getitem__(self, i):
		return self._rows[i]

	def get_all_rows(self):
		return [r for r in self._rows ]

	def get_rows(self, name='', exact=True):
		if name == '':
			exact=False
		if exact:
			return [r for r in self._rows if name == r.get_name() ]
		else:
			return [r for r in self._rows if re.search(name, r.get_name()) ]

class BEDrow:
	def __init__(self, s):
		self._contents = s.strip().strip(',').split('\t') 

	def __repr__(self):
		return ":".join(self._contents)
	
	def get_geneName(self):
		return self._contents[0]

	def get_name(self):
		return self._contents[1]

	def get_chr(self):
		return self._contents[2]

	def get_strand(self):
		return self._contents[3]

	def get_txStart(self):
		return int(self._contents[4])

	def get_txEnd(self):
		return int(self._contents[5])

	def get_txRange(self):
		return (self.get_txStart(), self.get_txEnd())

	def get_utr(self):
		result = [ (min(a,b), max(a,b)) for a, b in zip(self.get_txRange(), self.get_cdsRange()) ]
		for i in range(2):
			if result[i][0] == result[i][1]:
				result[i] = tuple()
		return result

	def get_cdsStart(self):
		return int(self._contents[6])

	def get_cdsEnd(self):
		return int(self._contents[7])

	def get_cdsRange(self):
		return (self.get_cdsStart(), self.get_cdsEnd())

	def get_exon_count(self):
		return int(self._contents[8])

	def get_exonStarts(self):
		return [int(a) for a in self._contents[9].split(',') if len(a) > 0]

	def get_exonEnds(self):
		return [int(a) for a in self._contents[10].split(',') if len(a) > 0]

	def get_exons(self):
		return [ (a,b) for a,b in zip( self.get_exonStarts(), self.get_exonEnds() ) ]


if __name__ == "__main__":

	bed = BED('../testing/data/refFlat.txt.gz.1')

	for i in range(5):
		r = bed[i]
		print r.get_name(), r.get_chr(), r.get_strand(), r.get_txRange(), r.get_cdsRange()
		print '\texons: ', r.get_exons() 
		print '\tUTR:   ', r.get_utr()
		print '\n'

	print bed.get_rows('NR_026820')
	
