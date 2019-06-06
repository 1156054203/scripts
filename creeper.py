#!/usr/bin/python
# -*- coding: utf-8 -*-
#get information for encode web

import sys
import time
import requests
from lxml import etree
from random import randint
import xlrd,xlwt

header = {
 'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.108 Safari/537.36',
 'Cookie':'gr_user_id=1f9ea7ea-462a-4a6f-9d55-156631fc6d45; bid=vPYpmmD30-k; ll="118282"; ue="codin; __utmz=30149280.1499577720.27.14.utmcsr=douban.com|utmccn=(referral)|utmcmd=referral|utmcct=/doulist/240962/; __utmv=30149280.3049; _vwo_uuid_v2=F04099A9dd; viewed="27607246_26356432"; ap=1; ps=y; push_noty_num=0; push_doumail_num=0; dbcl2="30496987:gZxPfTZW4y0"; ck=13ey; _pk_ref.100001.8cb4=%5B%22%22%2C%22%22%2C1515153574%2C%22https%3A%2F%2Fbook.douban.com%2Fmine%22%5D; __utma=30149280.833870293.1473539740.1514800523.1515153574.50; __utmc=30149280; _pk_id.100001.8cb4=255d8377ad92c57e.1473520329.20.1515153606.1514628010.'
}

base_url = "https://www.encodeproject.org/"
experiment = "%sexperiments/" % base_url
antibody = "%santibodies/" % base_url

workbook = xlrd.open_workbook(sys.argv[1])
sheet = workbook.sheet_by_index(1)

for cell in sheet.col_values(0):
    num = cell.strip()
    exper_url = experiment + num
    r1 = requests.get(exper_url, headers = header)
    tree1 = etree.HTML(r1.text)
    antibodies = tree1.xpath("//table[@class='table table-sortable']/tbody/tr[1]/td[5]/a/text()")
    control=tree1.xpath("//li[@class='multi-comma']/a/text()")
    if len(antibodies) == 0:
        sys.stdout.write("Error while get productid of %s\n")
        continue
    time.sleep(randint(1,3))
    anti = antibodies[0]
    anti_url = antibody + "/" + anti
    r2 = requests.get(anti_url, headers = header)
    tree2 = etree.HTML(r2.text)
    pro_id = tree2.xpath("//div[@data-test='productid']/dd/a/text()")
    sys.stdout.write("%s\t%s\t%s\t%s\n" % (num, control[0], anti, pro_id[0]))
    time.sleep(randint(1,3))
