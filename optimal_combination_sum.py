#!/usr/bin/python
#encoding:utf-8
from itertools import product
from itertools import combinations 
import itertools

num_list=[177,216,111,487,309,303,264,1245,228,220,106,31]
budget=1475
#min(summ-budget)


def test_combine(num_list):
    tmp_list=[]
    res_list=[]
    for i in range(len(num_list)+1):
        tmp_list+=list(combinations(num_list, i))
    for x in tmp_list:
        a=sum(x)-budget
        if  a > 0:
            res_list.append(list(x))
    n=len(res_list)
    mi=res_list[0]
    for a in range(n-1):
        if sum(res_list[a]) < sum(mi):
            mi=res_list[a]
    print (mi,sum(mi)) 


test_combine(num_list)
