# -*- encoding: utf-8 -
# env: python3.7

import os
import time
import json
import requests
import sqlite3
from pprint import pprint
# from pyquery import PyQuery as pq


def ncbi_search(keywords:dict={}, do_type='find/detail/summary'):
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    search_url = base_url + 'esearch.fcgi'      # 查询
    detail_url = base_url + 'efetch.fcgi'       # 详细信息
    summary_url = base_url + 'esummary.fcgi'     # 查询简要信息
    # keywords = {
    #     'db': 'gene',        # 数据库
    #     'term': 'CUL7+mouse',      # 查询关键词
    #     'RetMax': '100',     # 单页最大数量
    #     'id': '7157',        # 输入的ID
    #     'sort': 'date',      # 排序规则
    # }
    keywords['api_key'] = 'f2f08fd5065e0234025a65ddf48fdfd6a408'   # NCBI api_key 139563281@11.com
    if do_type == 'find':
        url = base_url + 'esearch.fcgi'
    elif do_type == 'detail':
        url = base_url + 'efetch.fcgi'     # 查询简要信息
    elif do_type == 'summary':
        url = base_url + 'esummary.fcgi'       # 详细信息'
    result = requests.get(url, params=keywords)
    # if result:
    #     result = result.json()
    # else:
    #     result = '无返回结果，参数错误'
    
    return result.text


def main():
    # # 搜索文献
    # keywords = {
    #     'db': 'pubmed',        # 数据库
    #     'term': 'CUL7+mouse',      # 查询关键词
    #     'RetMax': '100',     # 单页最大数量
    #     'sort': 'date',      # 排序规则
    # }
    # result_articel = ncbi_search(keywords, do_type='find')
    # # 获取文献ID
    # ids = result_articel['esearchresult']['idlist']
    # # 获取文献详细信息
    # result_detail = ncbi_search(ids, do_type='detail')
    # print('\n=====detail=====\n')
    # pprint(result_detail)
    # 获取文献简要信息
    # result_summary = ncbi_search(ids, do_type='summary')
    # print('\n=====result_summary=====\n')
    # pprint(result_summary)


    # # 搜索基因
    # keywords = {
    #     'db': 'gene',        # 数据库
    #     'term': 'CUL7',      # 查询关键词
    #     'RetMax': '100',     # 单页最大数量
    #     'sort': 'date',      # 排序规则
    # }
    # result_gene_find = ncbi_search(keywords, do_type='find')
    # print('\n=====gene_find=====\n')
    # pprint(result_gene_find)


    # 获取某基因的简要信息
    keywords = {
        'db': 'gene',        # 数据库
        'term': 'CUL7',      # 查询关键词
        'RetMax': '100',     # 单页最大数量
        'id': 7157,
        'sort': 'date',      # 排序规则
    }
    result_gene = ncbi_search(keywords, do_type='detail')
    print('\n=====gene_detail=====\n')
    # pprint(result_gene)
    with open('result_gene.json', 'w', encoding='utf8') as fp:
        fp.write(result_gene)

    # # 某文章的详细信息
    # keywords = {
    #     'db': 'pubmed',        # 数据库
    #     'term': 'CUL7',      # 查询关键词
    #     'retmode': 'xml',
    #     'rettype': 'xml',     
    #     'id': '35305671,33618520',
    #     'sort': 'date',      # 排序规则
    # }
    
    # result_article = ncbi_search(keywords, do_type='detail')
    # print('\n=====article_detail=====\n')
    # pprint(result_article)


if __name__ == '__main__':
    main()