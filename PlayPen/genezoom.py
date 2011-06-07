'''
Created on Jun 7, 2011

@author: jcc7
'''
import vcf

class CrossTable:
    def __init__(self, a, b):
        self.akeys=list(set(a))
        self.bkeys=list(set(b))
        self.akeys.sort()
        self.bkeys.sort()
        self.table = [[[] for bi in range(len(self.bkeys))] for ai in range(len(self.akeys))]
        for i in range(len(a)):
            indexa=self.akeys.index(a[i])
            indexb=self.bkeys.index(b[i])
            try:
                self.table[indexa][indexb]=self.table[indexa][indexb]+1
            except TypeError:
                self.table[indexa][indexb]=1
    def getakeys(self):
        return self.akeys
    def getbkeys(self):
        return self.bkeys
    def valueAt(self, i, j):
        return self.table[i][j]

#receives a pair of vectors and creates a table of data from them, based on a list of dictionaries
def xtabs(a, b):
    #create a list of dictionaries
    data=[{}]
    #go through each entry in vector a
    for i in a:
        nonexistant=True
        #go through each entry in data, checking if the entry in vector a exists in our data
        for entry in data:
            if i==entry:
                nonexistant=False
                if data[entry].has_key(b[i]):
                    temp=data[entry][b[i]]
                    data[entry][b[i]]=temp+1
                else:
                    data[entry][b[i]]=1
        #if this value does not yet exist in data, create a dictionary with a proper entry and append it to our list
        if nonexistant:
            temp={}
            temp[b[i]]=1
            data.append(temp)
    return data


#testing xtabs
if __name__ == '__main__':
    vecA=('case', 'case', 'control', 'control', 'case', 'control', 'case', 'case', 'case')
    vecB=('2', '1', '0', '0', '0', '1', '2', '1', '2')
    values=CrossTable(vecA, vecB)
    for i in range(2):
        for j in range(3):
            print values.valueAt(i, j)