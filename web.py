#!/usr/bin/env python3
#coding:utf-8
import os
import re
import sys
from bs4 import BeautifulSoup
from selenium import webdriver
if os.path.exists('/path'):
    sys.path.append('/path')
else:
    sys.path.append('.')
import common

def get_html(*argv_infos):
    '''
    web.get_html
    :param argv_infos: url+...
    :return: html_file
    '''
    # 1 set
    argv_dict = common.argv(*argv_infos)
    if 'url' not in argv_dict:
        print('Url Not Known!'); print(argv_infos)
        sys.exit()

    #2
    web_driver = webdriver.Chrome()
    web_driver.get(argv_dict['url'])
    html_page = BeautifulSoup(web_driver.page_source)
    print(html_page.prettify())

if __name__ == '__main__':
    if re.search('^(get_html)$',sys.argv[1],re.IGNORECASE):
        print( get_html(*sys.argv[2:]) )
    else:
        print('Ex:')
        print('\tpython3 *py get_html url+...;')
    sys.exit()








