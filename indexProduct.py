#!/env/python
# -*- coding: UTF-8 -*-

#输入index文件,输出与其相距为n的全部序列
from itertools import product
import sys

##计算两个字符串的汉明距离
def hanmingDistance(s1,s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return (sum(x1 != x2 for x1,x2 in zip(s1,s2)))

def diff_(s,n):
    s = s.upper()
    items = ["A","G","C","T"]

    lst_ = []
    for p in product(items,repeat=len(s)):
        lst_.append("".join(p))
        
    lst2_ = []
    for seq in lst_:
        if hanmingDistance(s,seq) == n:
            lst2_.append(seq)
    return lst2_
s = "AGCT"

idx_lst = []
with open(sys.argv[1],"r") as fi:
    for line in fi:
        line = line.strip("\n")
        if not line.startswith("#"):
            comb_idx = line.split("=")[0]+line.split("=")[1]
            idx_lst.append(comb_idx)
result = []
for idx in idx_lst:
    result.append(diff_(idx[-6:],1))

print(result)
