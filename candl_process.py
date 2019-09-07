#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-07 11:15:12
Last Modified time: 2019-09-07 11:15:12
Description: Extract non-synonymous mutations in the CanDL dataset. Since the genome information in this file have wrong info, so we get it by transvar tool.
'''

import pandas as pd
import os

path = os.getcwd()+"/output/"
data = "candl"
chrcol,poscol,refcol,altcol,aarefcol,aaaltcol,aaposcol = 'chr','start','ref','alt','Normal AA','Mutation AA','Peptide Position'
genecol,aacol = 'gene','aachange'
base = ['gene','chr','start','ref','alt','aaref','aaalt','aapos','aachange','mutationtype','key','index','mutation_ref']

def add_mutatype_drop_col():
	candl = pd.read_csv("/data1/wangs/PremPTS/data/CanDL/candl-results-20190709231125.csv",dtype=str)
	# delete no amino acid mutations
	candl = candl[candl['Normal AA'].notnull()] 
	candl['aachange'] = candl['Normal AA']+candl['Peptide Position']+candl['Mutation AA']
	candl['aachange'].replace('X','*', inplace=True)
	candl['key']= candl['Gene']+':'+candl['aachange']
	candl.loc[candl['Mutation AA']=='X','mutationtype'] = 'nonsense'
	candl.loc[(candl['Mutation AA']!='X')&(candl['Normal AA']!=candl['Mutation AA']),'mutationtype'] = 'missense' 
	candl.to_csv(path + "candl_add_mutatype.tsv", sep='\t',index=False)
	## genome of this file have wrong info, so we get it by transvar tool
	candl.drop(['ID','Chromosome','DNA Position','Codon','Transcript','Gene Strand','RNA Position','Exon','Mutation Codon'],axis=1,inplace=True)
	candl = candl.drop_duplicates()
	candl.rename(columns={'Gene':'gene','Normal AA':'aaref','Mutation AA':'aaalt','Peptide Position':'aapos'},inplace=True)
	candl.to_csv(path + "candl_del_cols.tsv", sep='\t',index=False)


def transvar_input():
	f = pd.read_csv(path + "candl_del_cols.tsv", sep='\t',dtype=str)
	f['key'].drop_duplicates().to_csv(path+"candl_transvar_input.txt", sep='\t',index=False)


def get_transvar():
	os.system('transvar panno -l /data1/wangs/PremPTS/data/CanDL/output/candl_transvar_input.txt --ensembl --refversion hg19 > /data1/wangs/PremPTS/data/CanDL/output/candl_transvar_output.txt 2> transvar.error')


def process_transvar():
	f = pd.read_csv(path + data + '_transvar_output.txt',sep='\t',dtype=str)
	# extract non single base mutation and can transaction by transvar
	f = f[(f['transcript'] != '.') & (~f['coordinates(gDNA/cDNA/protein)'].str.contains('del')) & (~f['info'].str.contains('Unclassified'))]
	# we only select the first output of transvar results
	f = f.drop_duplicates(['input'])
	f['gDNA'] = f['coordinates(gDNA/cDNA/protein)'].str.split('/').str[0]
	f['Mutation CDS'] = f['coordinates(gDNA/cDNA/protein)'].str.split('/').str[1]
	f['chr'] = f['gDNA'].str.split(':').str[0].str.strip('chr')
	f['start'] = f['gDNA'].str.findall('\d{3,}').str[0]
	f['ref'] = f['gDNA'].str[-3]
	f['alt'] = f['gDNA'].str[-1]
	f.to_csv(path + data + '_transvar_add_info.tsv', sep='\t', index=False)


def add_transvar_info():
	f1 = pd.read_csv(path + 'candl_del_cols.tsv', sep='\t',dtype=str)
	# add transvar chrmosome information
	transvar = pd.read_csv(path + data + '_transvar_add_info.tsv', dtype=str, sep='\t')
	transvar['gene'] = transvar['input'].str.split(':').str[0]
	transvar['aachange'] = transvar['input'].str.split(':').str[1]
	transvar.drop(['transcript'], axis=1, inplace=True)
	fw = f1.merge(transvar, how='inner', on=['gene', 'aachange'])
	fw['index'] = fw['chr'] + ':' + fw['start'] + ':' + fw['ref'] + ':' + fw['alt'] + ':' + fw['aaref'] + ':' + fw['aaalt']
	fw['mutation_ref'] = fw['gene']+';chr'+fw['chr'] + ';' + fw['start'] + ';'+fw['start'] +';'+ fw['ref'] + '>' + fw['alt']
	fw.to_csv(path + data + '_merge_transvar_info.tsv', sep='\t', index=False)


def clean():
    f = pd.read_csv(path + data + '_merge_transvar_info.tsv', sep='\t', dtype=str)
    need = f[base+['Level of Evidence']].drop_duplicates()
    need['source'] = 'CanDL'
    need.to_csv(path + "candl_ns_clean.tsv",sep='\t',index=False)


def main():
	add_mutatype_drop_col()
	transvar_input()
	get_transvar()
	process_transvar()
	add_transvar_info()
	clean()


if __name__ == '__main__':
	main()
