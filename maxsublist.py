#!/bin/env python

class maxsum:
    def __init__(self,list):
        self.list=list
    def getmaxsub(self):
        max= 0
        x,y=0,0
        sublist=[]
        for i in range(0,len(self.list)+1):
         for j in range(0,len(self.list)+1):
          res = sum(self.list[i:j])
          if res > max:
           max = res
           x = i
           y = j
           sublist=self.list[i:j]
        return max,sublist

    def sublist(self):
        sum = self.list[0]
        presum = 0
        for i in self.list:
            if presum < 0:
                presum = i
            else:
                presum += i
            sum = max(presum,sum)
        return sum
if __name__=='__main__':
    list=[-2, 1, 7, -4, 5, 2, -3, -6, 4, 3, -8, -1, 6, -7, -9, -5]
    a=maxsum(list)
    print('max: %s' % a.sublist())
    print('max: %s' % a.getmaxsub()[0],'sublist: %s' % a.getmaxsub()[1])
