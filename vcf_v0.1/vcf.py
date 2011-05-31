###
### Version 0.1b
### vcf.py - a prototype Python API for VCF parsing
###
### Written by Mikhail Spivakov 2011, current version distributed under the terms of CC BY-SA licence (Creative Commons Attribution-ShareAlike)
### Please report bugs/comments/suggestions to spivakov@ebi.ac.uk
### 
###

import re
import sys
import cPickle
import random
import pprint

# Creates subdictionaries from dictionaries containing only a specified subset of the original keys  
extract = lambda keys, dict: reduce(lambda x, y: x.update({y[0]:y[1]}) or x, map(None, keys, map(dict.get, keys)), {})

class CoordData:
	data = {}
	sk = set()

	def readFromPickle (self, fname):
		print >> sys.stderr, 'Unpickling coordData...'
		inputfile = open (fname)
		self.data = cPickle.load (inputfile)		
		inputfile.close()
		self.sk = set (self.data.keys()) # so much faster than lists and tuples for membership checks!

	def pickle (self, fname):
		print >> sys.stderr, 'Pickling coordData...'
		outputfile = open (fname, "w")
		cPickle.dump (self.data, outputfile, True)
	 	outputfile.close()

	def get(self, coordslist): # template for this method in derived classes
 		print >> sys.stderr, 'Warning: coordData\' original get() method invoked. This always returns an empty list' 
		return([])

	def __getitem__ (self, key):
		if type(key) is list or type(key) is set or type(key) is tuple:
			return self.get(key)
		else:
			return self.data[key]

	def __iter__ (self):
		self.current = 0
		return self

	def next (self):
		if self.current >= len(self.data.keys()):
			raise StopIteration
		else:
			self.current += 1
			res = self.data[self.data.keys()[self.current - 1]]
			res['COORD'] = self.data.keys()[self.current - 1]
			return res

	def items (self):
		return self.data.items()

	def keys (self):
		return self.data.keys()

	def values (self):
		return self.data.values()

	def __and__(self, other):
		return self.sk & other.sk

	def __or__(self, other):
		return self.sk | other.sk

	def __sub__(self, other):
		return self.sk - other.sk

	def __len__(self):
		return len(self.sk)

	def __repr__(self):
		if len(self)<20:
			return pprint.pformat(self.data)
		else:
			print "==Showing first 20 elements=="
			return pprint.pformat(dict((k, self.data[k]) for k in self.data.keys()[0:20])) 

	def __init__ (self):
		self.data = {}
		self.sk = set()

	def __del__ (self):
		self.data = {}
		self.sk = set()


def getStrain (tupleOfDicts, sampleNo, field="NUM"):
# given a tuple of dictionaries that have a field (def="NUM")
# extract the dictionary (assuming it's unique!) with dict[field] == sampleNo
        i=0
        ok=0
        for sample in tupleOfDicts:            
	        if sample['NUM'] == sampleNo:
        		ok = 1
                        break
                i+=1
        if ok:
        	return tupleOfDicts[i]
	else:
		return None

def coord (c):
        ccoord = re.search("(.+)@(.+)", c)
        return (ccoord.group(1), ccoord.group(2))

class VCF(CoordData):
        # note CoordData.data in this case is a dictionary of dictionaries 
	# containing pool-wide tags (including also "ALT" and "QUAL")) plus 
	# sample-specific info that is stored as a tuple of dictionaries (one tuple element per sample) under the '_SAMPLES' key
  
	# inherits CoordData operators:
	# -- coords level: y or z; y and z; y - z
	# -- slices: y[pos]; y[[pos1,pos2,...,posN]]
	# -- iterator:
	#     for x in y:
	# 	  x.keys()[0] # pos
	#	  x.values[0] # all the rest
	
	# Important new operators and differences from CoordData and its other derivatives:
	#
	# -- filter:
	# 	y.filter("COORD in set(['chrX@100', 'chr3@23234'])") or,if the condition includes at least one occurence of '>', '<', or '=', can use a shorthand:
	#	y["QUAL > 50 and DP < 1 and any('GQ>2')"]
	#    -- sample-specific tags are accessible via:
	#		any(cond), all(cond), atleast(N, cond), forsample(sampleNUM, cond)
	#    -- to compare between sample-specific tags:
	#		crosscomp(cond), where each sample-specific tag becomes a dictionary with sampleNo as key !
	#		  Example: y.filter("crosscomp('GT1[1]==GT1[2]')")
	# Remember that samples have an aux tag, NUM, that tells which ones they are in the order of columns in VCF 
	#
	# -- The result of get(), filter() and [] is a class!
	#	this is so we can do y.filter["AD>12"].field("DP")
	#
	# -- get field:
	#	y.field 

	def readFromFile(self, fname, refcol=3, altcol=4, qualcol=5, infocol=7, sampleformatcol=8, checkNonRef = False, commentChar = '#', verbose = True, readFromList = False):
		      # all col numbers are 0-based
			self.data = {}
			if verbose:
				print >> sys.stderr, 'Reading variation data...'
			if readFromList:
				snpfile = fname
			else:
				snpfile = open (fname)
			i=0
			for line in snpfile:
				line = line.strip()
				if line[0] == '#':
					continue
				linearr = line.split("\t")
				if checkNonRef:
					if linearr[refcol].upper() == linearr[altcol].upper():
						continue
				pos = linearr[0] + "@" + linearr[1]

				info = {'ALT' : linearr[altcol], 'REF' : linearr[refcol], 'QUAL' : float(linearr[qualcol])}

				infoarr = linearr[infocol].split(";")
				for field in infoarr:
					ft = field.split("=")
					tag = ft[0]	
					if re.search("\.", ft[1]):
						info[tag] = float(ft[1])
					else:
						info[tag] = int(ft[1])
				
											
					samlist = []
					strfields = linearr[sampleformatcol].split(":")
					samcount = 1
					for col in range(sampleformatcol+1, len(linearr)):
						strarr = linearr[col].split(":")
		
						if strarr[0] == "./.":
							samcount += 1
							continue
		
						thissam = {'NUM' : samcount}
						tagi = 0
						for tag in strfields:
							#print strarr, tagi, tag, strarr[tagi]
							if tag=="GT":
								gts = re.split("\||\/", strarr[tagi])				
								for gi in range(0, len(gts)):
									if gts[gi]=="0":
										thissam['GT%d' % (gi+1)] = "ref"
									elif gts[gi]=="1":
										thissam['GT%d' % (gi+1)] = "alt"
							elif tag=="AD":
								ads = strarr[tagi].split(",")
								thissam['ADref'] = int(ads[0])
								thissam['ADalt'] = int(ads[1])
							elif tag=="GL":		
								gls = strarr[tagi].split(",")
								thissam['GLrr'] = float(gls[0])
								thissam['GLra'] = float(gls[1])
								thissam['GLaa'] = float(gls[2])
							else:
								if re.search("\.", strarr[tagi]):
								    thissam[tag] = float(strarr[tagi])
								else:
								    thissam[tag] = int(strarr[tagi])
							tagi = tagi+1
						samlist.append(thissam)
						samcount += 1
					info['_SAMPLES'] = tuple(samlist)

				self.data[pos] = info

				i = i+1
				if (not (i+1) % 100000) and verbose:
					print >> sys.stderr, '.',

			self.sk = set (self.data.keys())
			if not readFromList:
				snpfile.close()
			if verbose:
				print >> sys.stderr, '\nDone!\n'

	def readFromList(self, list, verbose = False, *vcfFormatArgs):
			self.readFromFile(list, readFromList=True, verbose = verbose, *vcfFormatArgs)

	def get(self, coordslist, verbose = False):
	# Important! Returns a VCF class object and not a dictionary!

			keysInCoordsList = [thiskey for thiskey in coordslist if thiskey in self.sk]
			vcfInCoordsList = extract(keysInCoordsList, self.data) # note that the result is a dictionary of tuples - just like self.data!

			if verbose:
				print >> sys.stderr, 'len keysInCoordsList = ', len(keysInCoordsList)

			res = VCF()
			res.data = vcfInCoordsList
			res.sk = set(vcfInCoordsList.keys())
			return res

	def filter(self, expr):
	# expressions can use COORD, CHR, POS, REF, ALT and QUAL as well as any tag from the INFO field (in upper case!)
	# additionally, for the level of samples, it can be any("condition"), all("condition"), atleast(N, "condition") or forsample(sampleN, condition)
	# note samples with "./." for a SNP are not included (it may cause more confusion if they are...)

			pos = 0 # stub; should be defined below when iterating through self
			this = {} # same as above

			def __populatens(sample):
				ns = {}
				for (tag, val) in sample.items():
					ns[tag] = val
				return ns

			def any(condition):
				# note pos and this come from the scope of self.filter()
				for sample in this['_SAMPLES']:
					ns = __populatens(sample)
					if eval(condition, ns):
						return True
				return False
						
			def all(condition):
				# note pos and this come from the scope of self.filter
				for sample in this['_SAMPLES']:
					ns = __populatens(sample)
					if not eval(condition, ns):
						return False
				return True

			def atleast(n, condition):
				passed = 0
				for sample in this['_SAMPLES']:
					ns = __populatens(sample)
					if eval(condition, ns):
						passed += 1
				return (passed >= n)

			def forsample(sampleNo, condition):
				# will return false if sampleNo does not exist
				sample = getStrain(this['_SAMPLES'], sampleNo)
				if not sample:
					return False
				ns = __populatens(sample)
				return eval(condition, ns) 

			def crosscomp(condition):
				ns = {}
				for sample in this['_SAMPLES']:
					num = sample['NUM']
					for (tag, val) in sample.items():
						if tag=='NUM':
							continue
						if not tag in ns:
							ns[tag] = {}
						ns[tag][num] = val
				try:
					return eval(condition, ns)
				except (KeyError, NameError):
					return False

			which = []
			for (pos, this) in self.items():
				ns = {'any' : any, 'all' : all, 'atleast' : atleast, 'forsample' : forsample, 'crosscomp' : crosscomp}
				ns['COORD'] = pos
				ccoord = re.search("(.+)@(.+)", pos)
				ns['CHR'] = ccoord.group(1)
				ns['POS'] = ccoord.group(2)
				for (tag, val) in this.items():
					ns[tag] = val
				if eval(expr, ns):
					which.append(ns['COORD'])
			return self.get(set(which))


	def field(self, field, sampleNo = None):
			res = []
			if field=="COORD":
				return self.keys()
			elif field=="CHR":
				return [re.search("(.+)@", coord).group(1) for coord in self.keys()]
			elif field=="POS":
				return [re.search("@(.+)", coord).group(1) for coord in self.keys()]
			elif not sampleNo:
				for this in self.values():
					res.append(this[field])
			else:
				for this in self.values():
					# returns -1 if sample not found
					thisStrain = getStrain(this['_SAMPLES'], sampleNo)
					if thisStrain:
						res.append(thisStrain[field])
			return res

	def __getitem__ (self, key):
			if type(key) is list or type(key) is set or type(key) is tuple:
				return self.get(key)
			elif (type(key) is str and re.search("[><=]", key)):
				return self.filter(key)
			else:
				res = VCF()
				res.data = {key: self.data[key]}
				res.sk = set([key,])
				return res

def VCFfilter(fname, filter=None, field=None, count=False, chunk = 10000, maxchunks = None, *vcfFormatArgs):
	# Reads a VCF file chunk by chunk, performing filtering and field extraction (or element counting) for each chunk.
	#
	# Output:
	# if field==None and count==False: a VCF object containing all elements passing the filter (or all elements if filter is not defined),
	# if field is specified: a list containing this field for all elements passing the filter,
	# if count is specified: an integer corresponding to the total number of elements passing the filter.
	#
	# With a test VCF consisting of 50k lines, the optimum chunk size was ~10000. 
	# The limiting stage seems to be the merger of dictionaries, so it may actually perform faster with filtering than without it.

		if field and count:
			raise ValueError('When count is True, field should be None, because only the count of entries passing the filter will be returned')
		f = open(fname)
		filtered = {} # filtered VCF.data
		yff = [] # (filtered) VCF.data field
		counter = 0
		nchunks = 0
		while True:
			l=[]
			i=0
			stp=False
			while i<chunk:
				line = f.readline()
				if len(line):
					l.append(line)
					i+=1
				else:
					stp = True
					break

			print >> sys.stderr, ".",
			nchunks += 1
			if nchunks >= maxchunks and maxchunks:
				stp = True

			y = VCF()
			y.readFromList(l, *vcfFormatArgs)
			if filter:
				y = y.filter(filter)
                
			if count:
        	                counter+=len(y)
			elif field:
				if type(field)==str:
					yff.extend(y.field(field))
				elif len(field)==2:
					yff.extend(y.field(field[0], field[1]))
				else:
					raise ValueError('field can either be a string or a tuple of size 2: (string, sampleNo)')
			else:
                        	filtered = dict(filtered, **y.data)

			l = []
			if stp:
				break
		f.close()
		print >> sys.stderr, ""

		if count:
			return counter
		elif field:
			return yff
		else:
			z = VCF()
			z.data = filtered
			z.sk = set(z.data.keys())
			return (z)
