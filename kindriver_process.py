#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-07 11:16:24
Last Modified time: 2019-09-07 11:16:24
Description: Extract non-synonymous mutations in the KinDriver dataset. The information of genome change are get by TransVar tool.
'''

import os
import re
import pandas as pd

path = '/data1/wangs/PremPTS/data/KinDriver/'
outpath = '/data1/wangs/PremPTS/data/KinDriver/output/'
data = 'kindriver'
chrcol,poscol,refcol,altcol,aarefcol,aaaltcol,aaposcol = 'chr','start','ref','alt','aaref','aaalt','aapos'
genecol,aacol = 'gene','aachange'
base = ['gene','chr','start','ref','alt','aaref','aaalt','aapos','aachange','mutationtype','key','index','mutation_ref']


def uniprot_id():
	kin = pd.read_csv(path + 'kindriver_simple_mutations_v82.txt', sep='\t', dtype=str)
	kin['uniprot_id'].drop_duplicates().to_csv(path + "kindriver_v82_uniprotID.txt", sep='\t', index=False)


def uniprot_to_gene_ns():
	kin = pd.read_csv(path + 'kindriver_simple_mutations_v82.txt', sep='\t', dtype=str)
	map = pd.read_csv(outpath + 'kindriver_v82_uniprotID_to_gene.txt', sep='\t', dtype=str)
	kin = kin.merge(map, how='left', left_on=['uniprot_id'], right_on=['From']).drop(['From'],axis=1)
	kin.rename(columns={'To':'gene','mutation':'aachange','mut_type':'mutationtype'},inplace=True)
	# non-synonymous
	ns = kin[kin['mutationtype'].isin(['missense','nonsense'])]
	# 一行中有多个突变的分开展示
	ns = ns.drop('aachange', axis=1).join(ns['aachange'].str.split('+', expand=True).stack().reset_index(level=1, drop=True).rename('aachange'))
	ns['aaref'] = ns['aachange'].str.split('\d+').str[0]
	ns['aaalt'] = ns['aachange'].str.split('\d+').str[1]
	ns['aaalt'].replace('*','X',inplace=True)
	ns['aapos'] = ns['aachange'].str.findall('\d+').str[0]
	ns['key'] = ns['gene']+':'+ns['aachange']
	ns.to_csv(outpath + "kindriver_simple_mutations_ns_v82.tsv",sep='\t',index=False)


def transvar_input():
	f = pd.read_csv(outpath + "kindriver_simple_mutations_ns_v82.tsv",sep='\t',dtype=str)
	f['key'].drop_duplicates().to_csv(outpath+"kindriver_transvar_input.txt", sep='\t',index=False)


def get_transvar():
	os.system("transvar panno -l /data1/wangs/PremPTS/data/KinDriver/output/kindriver_transvar_input.txt --ensembl --refversion hg19 > /data1/wangs/PremPTS/data/KinDriver/output/kindriver_transvar_output.txt 2> kindriver.error")


def process_transvar():
	f = pd.read_csv(outpath + data + '_transvar_output.txt',sep='\t',dtype=str)
	# extract non single base mutation and can transaction by transvar
	f = f[(f['transcript'] != '.') & (~f['coordinates(gDNA/cDNA/protein)'].str.contains('del')) & (~f['info'].str.contains('Unclassified'))]
	f = f.drop_duplicates(['input'])
	f['gDNA'] = f['coordinates(gDNA/cDNA/protein)'].str.split('/').str[0]
	f['Mutation CDS'] = f['coordinates(gDNA/cDNA/protein)'].str.split('/').str[1]
	f['chr'] = f['gDNA'].str.split(':').str[0].str.strip('chr')
	f['start'] = f['gDNA'].str.findall('\d{3,}').str[0]
	f['ref'] = f['gDNA'].str[-3]
	f['alt'] = f['gDNA'].str[-1]
	f.to_csv(outpath + data + '_transvar_add_info.tsv', sep='\t', index=False)


def add_transvar_info():
	f1 = pd.read_csv(outpath + 'kindriver_simple_mutations_ns_v82.tsv',sep='\t',dtype=str)
	# add transvar chrmosome information
	transvar = pd.read_csv(outpath + data + '_transvar_add_info.tsv', dtype=str, sep='\t')
	transvar[genecol] = transvar['input'].str.split(':').str[0]
	transvar[aacol] = transvar['input'].str.split(':').str[1]
	fw = f1.merge(transvar, how='inner', on=[genecol, aacol])
	fw = fw.drop(['transcript'], axis=1)
	fw['index'] = fw['chr'] + ':' + fw['start'] + ':' + fw['ref'] + ':' + fw['alt'] + ':' + fw['aaref'] + ':' + fw['aaalt']
	fw['mutation_ref'] = fw['gene']+';chr'+fw['chr'] + ';' + fw['start'] + ';'+fw['start'] +';'+ fw['ref'] + '>' + fw['alt']
	fw.to_csv(outpath + data + '_merge_transvar_info.tsv', sep='\t', index=False)


def clean():
    f = pd.read_csv(outpath + data + '_merge_transvar_info.tsv', sep='\t', dtype=str)
    need = f[base+['validation']].drop_duplicates()
    need['source'] = 'KinDriver'
    need.to_csv(outpath + data + "_ns_clean.tsv",sep='\t',index=False)


def get_score():
	os.chdir('/data/wangs/MutaAnalysis')
	os.system("python data_process.py -p /data/wangs/MutaAnalysis/Hpm/KinDriver/output/ -d kindriver -i kindriver_ns_clean.tsv")


if __name__ == '__main__':
	uniprot_id()
	uniprot_to_gene_ns()
	transvar_input()
	get_transvar()
	process_transvar()
	add_transvar_info()
	clean()