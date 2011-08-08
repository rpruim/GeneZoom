import re
#import os
#import sys
import string
import logging


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
		return int(self._items[1])

	def get_locus(self):
		return ( int(self._items[0]),  int(self._items[1]) )

	def get_name(self):
		return self._items[2]

	def get_refAllele(self):
		return self._items[3]

	def get_altAllele(self, key=None):
		return getItemByNameFromEqList( self._items[4].split(',') , key )

	def get_numAlleles(self):
		return 1 + len(self.get_altAllele())

	def is_indel(self):
		logging.debug('Could implement a better check for indels')
		return not len(self.get_refAllele()) == len(self.get_altAllele()[0])

	def get_qual(self, key=None):
		return getItemByNameFromEqList( self._items[5].split(';') , key )

	def get_filter(self, key=None):
		return getItemByNameFromEqList( self._items[6].split(';') , key )

	def get_info(self, key=None):
		return getItemByNameFromEqList( self._items[7].split(';') , key )

	def get_format(self, key=None):
		return getItemByNameFromEqList( self._items[8].split(';') , key )
	'''
	From VCF 4.0 spec:
	If genotype information is present, then the same types of data must be present
	for all samples. First a FORMAT field is given specifying the data types and
	order. This is followed by one field per sample, with the colon-separated data
	in this field corresponding to the types specified in the format. The first
	sub-field must always be the genotype (GT).
	'''
	def get_genotypes(self):
		return self.get('GT') 
		# start=9
		# return [ g.split(":")[0] for g in self._items[start:] ]

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

	def get(self, key, cast=str):
		i = min( [ a for (a,b) in zip( range(len(self._items[8])), self._items[8].split(':')) if b == key ] )
		result = []
		for item in self._items[9:]:
			try:
				result.append( cast(item.split(':')[i]) )
			except IndexError:
				result.append(None)
		return result

	def __len__(self):
		return len(self._items)

	def __repr__(self):
		return string.join(self._items, "\t")

	def checkFilter(self, userFilter):
		for f in self.get_filter():
			if f=="PASS":
				return True
			elif userFilter!=None:
				for uf in userFilter:
					if uf == f:
						return False
			else:
				return False
		return True