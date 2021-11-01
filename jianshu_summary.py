#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys,argparse,requests
reload(sys)
sys.setdefaultencoding( "utf-8" )
from lxml import etree

def parse_args():
    parser = argparse.ArgumentParser(usage="""python %(prog)s -u url [-o output]""")
    parser.add_argument('-u',dest='url',metavar='url',help="jianshu link need to get summary information.")
    parser.add_argument('-o',dest='output',metavar='output',default='jianshu_summary.txt',help="the output file, default: %(default)s.")
    return parser

def get_article_summary(url,output):
    # 设置请求头
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.117 Safari/537.36'
    }
    r = requests.get(url + '?order_by=shared_at&page=1', headers=headers)
    dom = etree.HTML(r.text)

    #获取文章数量和最大页数
    article_num = int(dom.xpath('//div[@class="info"]//li[3]//p/text()')[0].strip())

    if (article_num % 9) == 0:
        max_page_num = article_num / 9
    else:
        max_page_num = (article_num / 9) + 1

    out = open(output, 'w')
    out.write('\t'.join(['link','title','read_num','comment_num','heart_num','diamond_num','release_time','\n']))
    xpath_items = '//ul[@class="note-list"]/li'
    for i in range(max_page_num):
        r = requests.get(url + '?order_by=shared_at&page=' + str(i + 1), headers=headers)
        article_list = etree.HTML(r.text).xpath(xpath_items)
        for article in article_list:
            link =  'https://www.jianshu.com' + ''.join(article.xpath('./div/a/@href')).strip()
            title = ''.join(article.xpath('./div/a/text()')).strip()
            read_num =  ''.join(article.xpath('.//div[@class="meta"]/a[1]/text()')).strip()
            comment_num = ''.join(article.xpath('.//div[@class="meta"]/a[2]/text()')).strip()
            heart_num = ''.join(article.xpath('.//div[@class="meta"]/span[2]/text()')).strip()
            diamond_num = ''.join(article.xpath('.//div[@class="meta"]/span[1]/text()')).strip()
            release_time = ''.join(article.xpath('.//div[@class="meta"]/span[3]/@data-shared-at')).strip()
            #print('link: %s, title: %s, read_num: %s, comment_num: %s, heart_num: %s, diamond_num: %s, release_time: %s' % (link,title,read_num,comment_num,heart_num,diamond_num,release_time))
            out.write('\t'.join([link,title,read_num,comment_num,heart_num,diamond_num,release_time,'\n']))
    out.close

if __name__ == "__main__":
    parse = parse_args()
    args = parse.parse_args()
    # Do not set Required =True, and print the complete help information without required parameters
    if not args.url:
        parse.print_help()
        sys.exit()
    get_article_summary(args.url,args.output)
