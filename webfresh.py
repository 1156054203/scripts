#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time,os,argparse
from random import randint
from selenium import webdriver
#chrome browser driver: http://chromedriver.storage.googleapis.com/index.html

def parse_args():
    usage = 'python webfresh.py -u url [-n num]'
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-u',dest='url',metavar='url',required=True,help="web link need to fresh.")
    parser.add_argument('-n',dest='num',metavar='num',type=int,default=200,help="fresh times, default: %(default)s.")
    return parser.parse_args()

def refresh(url,num):
    driver = webdriver.Chrome("/mnt/c/software/chrome/chromedriver.exe")
    driver.get(url)
    for i in range(num):
        time.sleep(randint(30,90))
        driver.refresh()
    driver.close()

if __name__ == "__main__":
    args = parse_args()
    refresh(args.url,args.num)
    print('fresh done....')
