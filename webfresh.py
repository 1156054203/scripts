#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time,os,argparse,sys
from random import randint
from selenium import webdriver
#chrome browser driver: http://chromedriver.storage.googleapis.com/index.html

def parse_args():
    parser = argparse.ArgumentParser(usage="""python %(prog)s -u url [-n num]""")
    parser.add_argument('-u',dest='url',metavar='url',help="web link need to fresh.")
    parser.add_argument('-n',dest='num',metavar='num',type=int,default=200,help="fresh times, default: %(default)s.")
    return parser

def refresh(url,num):
    driver = webdriver.Chrome("/mnt/c/software/chrome/chromedriver.exe")
    driver.get(url)
    for i in range(num):
        time.sleep(randint(30,90))
        driver.refresh()
    driver.close()

if __name__ == "__main__":
    parse = parse_args()
    args = parse.parse_args()
    # Do not set Required =True, and print the complete help information without required parameters
    if not args.url:
        parse.print_help()
        sys.exit()
    refresh(args.url,args.num)
    print('fresh done....')
