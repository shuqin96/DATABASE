#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-07 11:19:41
Last Modified time: 2019-09-07 11:19:41
Description: 
'''

import os
import pandas as pd

inpath = '/data/wangs/MutaAnalysis/Hpm/oncokb/'
outpath = '/data/wangs/MutaAnalysis/Hpm/oncokb/output/'

data = 'oncokb'
genecol = 'gene'
aacol = 'aachange'
base = ['gene','chr','start','ref','alt','aaref','aaalt','aapos','aachange','mutationtype','key','index','mutation_ref']

def process_oncokb():
	oncokb = pd.read_table(inpath + 'allAnnotatedVariantsOncoKB.txt')
	# single nucletiode mutation
	snv = oncokb[(oncokb['Alteration'].str.match('[A-Z]\d+[A-Z]$'))|(oncokb['Alteration'].str.match('[A-Z]\d+\*$'))]
	snv.rename(columns={'Protein Change':'aachange','Gene':'gene','Isoform':'transcript'},inplace=True)
	snv['key']=snv['gene']+':'+snv['aachange']
	snv.to_csv(outpath + 'oncokb_snv.tsv', sep='\t', index=False)
	snv['aapos'] = snv['aachange'].str.findall('\d+').str[0]
	snv['aaref'], snv['aaalt'] = snv['aachange'].str.split('\d+').str
	snv['aaalt'].replace('*', 'X', inplace=True)
	snv['mutationtype'] = ''
	snv.loc[snv['aaref'] == snv['aaalt'], 'mutationtype'] = 'silent'
	snv.loc[(snv['aaref'] != snv['aaalt']) & (snv['aaalt'] != '*'), 'mutationtype'] = 'missense'
	snv.loc[snv['aaalt'] == 'X', 'mutationtype'] = 'nonsense'
	ns = snv[snv['mutationtype'].isin(['missense','nonsense'])]
	ns.to_csv(outpath + 'oncokb_ns.tsv', sep='\t', index=False)


# Make transvar inputfile to get chr,ref,alt information
def transvar_input():
	f1 = pd.read_table(outpath + 'oncokb_ns.tsv')
	f1['key'].to_csv(outpath + 'oncokb_transvar_input.txt', sep='\t', index=False)


def get_transvar_out():
	os.chdir('/data/wangs/MutaAnalysis/Hpm/oncokb/output/')
	# 第一个>输出结果，第二个>输出报错信息
	os.system(('transvar panno -l %s --ensembl > %s 2> run_transvar.log') % ('oncokb_transvar_input.txt', 'oncokb_transvar_output.tsv')) 


def process_transvar():
	f = pd.read_table(outpath + data + '_transvar_output.tsv')
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
	f1 = pd.read_table(outpath + 'oncokb_ns.tsv')
	# add transvar chrmosome information
	transvar = pd.read_table(outpath + data + '_transvar_add_info.tsv', dtype=str)
	fw = f1.merge(transvar, how='inner', left_on=['key'], right_on=['input'])
	fw = fw.drop(['transcript','input'], axis=1)
	fw['index'] = fw['chr'] + ':' + fw['start'] + ':' + fw['ref'] + ':' + fw['alt'] + ':' + fw['aaref'] + ':' + fw['aaalt']
	fw['mutation_ref'] = fw['gene']+';chr'+fw['chr'] + ';' + fw['start'] + ';'+fw['start'] +';'+ fw['ref'] + '>' + fw['alt']
	fw.to_csv(outpath + data + '_merge_transvar_info.tsv', sep='\t', index=False)


def clean():
	f = pd.read_csv(outpath + data + '_merge_transvar_info.tsv', sep='\t', dtype=str)
	need = f[base+['Oncogenicity']].drop_duplicates()
	need['source'] = 'OncoKB'
	need.to_csv(outpath + data + "_ns_clean.tsv",sep='\t',index=False)


# we think mutations labeled as 'Oncogenic' or 'Likely Oncogenic' are more damage.
def damage():
    f = pd.read_csv(outpath +data +"_ns_clean.tsv",sep='\t',dtype=str)
    f = f[f['Oncogenicity'].isin(['Oncogenic','Likely Oncogenic'])]
    f.to_csv(outpath+data+"_ns_clean_damage_mutation.tsv",sep='\t',index=False)
   

def main():
	process_oncokb()
	transvar_input()
	get_transvar_out()
	process_transvar()
	add_transvar_info()
	clean()
	damage()


if __name__ == '__main__':
	main()
