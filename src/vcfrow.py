import re
#import os
#import sys
import string


def getItemByNameFromEqList(list, key):
	if key == None: 
		return list
	keys = [ a.split('=')[0] for a in list ]
	vals = [ a.split('=')[1] for a in list ]
	result = [ v for (k,v) in zip(keys, vals) if k == key ]
	return result

class VCFrow:
	def __init__(self, s):
		s.strip()
		self._items = s.split('\t')

	def get_chrom(self):
		return self._items[0]

	def get_pos(self):
		return self._items[1]

	def get_locus(self):
		return ( int(self._items[0]),  int(self._items[1]) )

	def get_name(self):
		return self._items[2]

	def get_refAllele(self):
		return self._items[3]

	def get_altAllele(self, key=None):
		return getItemByNameFromEqList( self._items[4].split(';') , key )

	def get_qual(self, key=None):
		return getItemByNameFromEqList( self._items[5].split(';') , key )

	def get_filter(self, key=None):
		return getItemByNameFromEqList( self._items[6].split(';') , key )

	def get_info(self, key=None):
		return getItemByNameFromEqList( self._items[7].split(';') , key )

	def get_format(self, key=None):
		return getItemByNameFromEqList( self._items[8].split(';') , key )

	def get_genotypes(self, start=9, stop=None):
		if stop == None:
			stop = len(self._items)
		return [ g.split(":")[0] for g in self._items[start:stop] ]
		# return [ string.split(g,":")[0] for g in self._items[start:stop] ]

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
		return string.join(self._items, "\t")
