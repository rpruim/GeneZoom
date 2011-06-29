#!/usr/bin/env python
'''
Created on Jun 7, 2011
Creates a cross table of two vectors, tabulating the number of 
occurrences of each pair x,y, with x a member of the first vector 
and y a member of the second
@author: jcc7
'''
#import vcf
#forms a table, of a list of lists
#case/control should be the first vector
#since it will be sorted, instances of two vectors of case/control and allele number will be formed in this format:
#        case     control
#table=[[ 23 ]    [  2  ]]
#       [ 15 ]    [  4  ]
#       [ 2  ]    [ 34  ]
class xTable:
    def __init__(self, a, b):
        self.akeys=list(set(a))
        self.bkeys=list(set(b))
        self.akeys.sort()
        self.bkeys.sort()
        self.table = [[0 for bi in range(len(self.bkeys))] for ai in range(len(self.akeys))]
        for i in range(len(a)):
            indexa=self.akeys.index(a[i])
            indexb=self.bkeys.index(b[i])
            self.table[indexa][indexb]=self.table[indexa][indexb]+1
    def getakeys(self):
        return self.akeys
    def getbkeys(self):
        return self.bkeys
    #returns a list containing the elements of the ith column
    def column(self, i):
        return self.table[i]
    #returns the value of the table at location ith column, jth row
    def valueAt(self, i, j):
        return self.table[i][j]

#testing xtabs
if __name__ == '__main__':
    #testing with junk data
    vecA=('case', 'case', 'control', 'control', 'case', 'control', 'case', 'case', 'case')
    vecB=('1/1', '0/1', '0/0', './.', './.', '1/0', '1/1', '1/0', '1/1')
    values=xTable(vecA, vecB)
    lista=values.getakeys()
    listb=values.getbkeys()
    #print out table
    print "\t",
    print listb,
    for i in range(len(lista)):
        print "\n"+lista[i]+"\t",
        for j in range(len(listb)):
            print "   %s" %values.valueAt(i, j),