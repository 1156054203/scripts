#!/usr/bin/python
# -*- coding: utf-8 -*-
#get information for NCBI GEO

import sys
import time
import requests
from lxml import etree

out=open('table.txt','w')

header = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/64.0.3282.186 Safari/537.36'}
url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100738'
response = requests.get(url,headers=header,params=None)
html=etree.HTML(response.text)
samples=html.xpath('//td[@valign="top"]/a/text()')[1:]
species=html.xpath('//tr[@valign="top"]/td/a/text()')[0]

table=[]
for sample in samples:
    req=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+sample,headers=header,params=None)
    html2=etree.HTML(req.text)
    infor=html2.xpath('//td[@style="text-align: justify"]/text()')[:4]
    title=infor[0]
    source=infor[1]
    celltype=infor[3].split(':')[1].strip()
    strain=infor[2].split(':')[1].strip()
    out.write(sample+'\t'+title+'\t'+source+'\t'+species+'\t'+celltype+'\t'+strain+'\n')

out.close()
