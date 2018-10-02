#!/env/python

##¿¿¿¿Lane¿index(6¿)¿¿¿¿

import itertools
import sys

##¿¿¿¿¿¿¿¿¿¿¿
def hanmingDistance(s1,s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return (sum(x1 != x2 for x1,x2 in zip(s1,s2)))
#print(hanmingDistance("ATCACG","CGATGT"))

##¿¿¿¿index¿¿¿¿¿¿¿¿
#sys.path.insert(0,r"C:\User\Desktop")
#idx_tuple = (0,)
idx_lst = []
with open(sys.argv[1],"r") as fi:
    for line in fi:
        line = line.strip("\n")
        if not line.startswith("#"):
            comb_idx = line.split("=")[0]+line.split("=")[1]
            idx_lst.append(comb_idx)

def comb(idx_lst):
    for c in itertools.combinations(idx_lst,2):
        c += (hanmingDistance(c[0][-6:],c[1][-6:]),)
        c = (c[0][:-6],c[0][-6:],c[1][:-6],c[1][-6:],c[2])
        print(c)

comb(idx_lst)
