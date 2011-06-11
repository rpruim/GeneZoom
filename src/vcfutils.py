###
### tools for dealing with vcf data
###

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
			result[ item ] = 0 

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
	filename = '../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz'

	v = tabix.Tabix(filename)

	if True:
		pos = 156705468 + 800000
		offset = 5000
		stuff = [row for row in v.query('1', pos, pos + offset)]
		print 'query(', pos, pos+offset, ')...'
		# print string.join( [str(l) for l in stuff], "\n" )
		for c in stuff:
			print  extract_chrom(c), extract_pos(c), ": ",
			print  tally( extract_genotypes(c) )
		print '\nTotal: ', len(stuff), 'markers'


