#!/usr/bin/env python
'''
Created on Jun 7, 2011
Creates a cross table of two vectors in the form of a dictionary of dictionaries, tabulating the number of 
occurrences of each pair x,y, with x a member of the first vector 
and y a member of the second.  Result is a two-dimensional dictionary with tabulation of occurrences of (x,y) in location table[x][y]
Also contains methods to ensure that both lists are the same size before making the cross table (methods must be called before creating the cross table).
'''
#import vcf
#forms a table, of a dictionary of dictionaries
#[case/control][allele count]

def cullList(master, keys, vals):
    '''keys and vals should be lists of the same length.  Returns a sorted and culled version of vals to match the keys in master rather than in keys
    In other words, it culls/expands and sorts vals, enabling a cross table to be made with a list of elements corresponding to master'''
    d=dict(zip(keys, vals)) 
    culledList=[d.get(s) for s in master]
    return culledList

class xTable:
    '''Creates a dictionary of dictionaries, tabulating occurrences of pairings of (x,y) in vectors a and b'''      
    def __init__(self, a, b):
        self.table = {}
        for (x, y) in zip(a, b):
            if self.table.has_key(x):
                if self.table[x].has_key(y):
                    self.table[x][y]=self.table[x][y]+1
                else:
                    self.table[x][y]=1
            else:
                self.table[x]={y:1}
                
  
    def printTable(self):
        '''Prints the table, in roughly table format'''
        print "\t",
        for entry in self.table.keys():
            print entry+"\t",
        print "\n",
        for entryb in self.table[entry]:
            print entryb+" \t ",
            for entryc in self.table.keys():
                if self.table[entryc].has_key(entryb):
                    print "%s\t "%self.table[entryc][entryb],
                else:
                    print "0\t ",
            print "\n",


    def getTable(self):
        '''Returns the table'''
        return self.table
    
#testing xtabs
if __name__ == '__main__':
    #testing with junk data
#    vecA=('case', 'case', 'control', 'control', 'case', 'control', 'case', 'case', 'case')
#    vecB=('1/1', '0/1', '0/0', './.', './.', '1/0', '1/1', '1/0', '1/1')
#    values=xTable(vecA, vecB)
#    values.printTable()
    master = ['1','2','3','4']
    keys = ['1','3','5','4','7']
    vals = ['one', 'three', 'five', 'four', 'seven']
    vals = ['odd', 'odd', 'odd', 'even', 'odd']
    print cullList(master, keys, vals)
    print xTable( master, cullList(master,keys,vals) ).getTable()
    print xTable( cullList(master,keys,vals), master ).getTable()
