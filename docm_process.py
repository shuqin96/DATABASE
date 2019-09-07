#!/usr/bin/python
# -*- coding: UTF-8 -*-


import os
import pandas as pd

inpath = "/data/wangs/MutaAnalysis/Hpm/docm3.2/"
outpath = '/data/wangs/MutaAnalysis/Hpm/docm3.2/output/'

base = ['gene','chr','start','ref','alt','aaref','aaalt','aapos','aachange','mutationtype','key','index','mutation_ref']

def a(row):  # output chr, ref, alt, aaref, aaalt
    return '{}:{}:{}:{}:{}:{}'.format(row[7], row[8], row[2], row[3], row[4], row[5])


def docmprocess1():  # output single nucleotide variants
	f1 = pd.read_table(inpath + "docm_version_3.2.tsv")
	f2 = f1[(f1['start'] == f1['stop']) & (f1['mutation_type'] != 'frameshift')]
	f2.to_csv(outpath + 'docm3.2_1.tsv', sep='\t', index=False)
	# Delete unuseful columns
	f3 = f2.drop(['pubmed_sources'], axis = 1)
	# rename
	f3.rename(columns={'read':'ref', 'variant':'alt','mutation_type':'mutationtype','chromosome':'chr'},inplace=True)
	f3.to_csv(outpath + 'docm3.2_2.tsv', sep='\t', index=False)


def docmprocess2():
	f1 = pd.read_table(outpath + 'docm3.2_2.tsv')
	silent = f1[f1['mutationtype'] == 'synonymous']
	silent.to_csv(outpath + 'docm3.2_Silent.tsv', sep='\t', index=False)
	non = f1[(f1['mutationtype'] == 'missense')|(f1['mutationtype'] == 'start_lost')|(f1['mutationtype'] == 'stop_lost')]
	non['mutationtype'].replace({'start_lost':'missense', 'stop_lost':'nonsense'}, inplace=True)
	non.to_csv(outpath + 'docm3.2_NS.tsv', sep='\t', index=False)


def docmprocess3():  # insert columns to write aa change, aaref, aaalt
    f1 = pd.read_table(outpath + "docm3.2_NS.tsv")
    f1['aachange'] = f1['amino_acid'].str.split('.').str[1]
    f1['aaref'] = f1['aachange'].str.split('\d+').str[0]
    f1['aaalt'] = f1['aachange'].str.split('\d+').str[1]
    f1['aapos'] = f1['aachange'].str.findall('\d+').str[0]
    f1['aaref'].replace('*', 'X', inplace=True)
    f1['aaalt'].replace('*', 'X', inplace=True)
    f1['key'] = f1['gene']+':'+f1['aachange']
    f1.to_csv(outpath + 'docm3.2_NS_insert_aa.tsv', sep='\t', index=False)


def docmprocess4():
	f1 = pd.read_table(outpath + "docm3.2_NS_insert_aa.tsv")
	acc = f1['hgvs'].str.split(':').str[0]
	f1.insert(9, 'transcript', acc)
	f1.to_csv(outpath + "docm3.2_NS_insert_accession.tsv", sep='\t', index=False)


def docmprocess5(): 
    f1 = pd.read_table(outpath + 'docm3.2_NS_insert_accession.tsv', dtype=str)
    f1['index'] = f1['chr'] + ':' + f1['start'] + ':' + f1['ref'] + ':' + f1['alt'] + ':' + f1['aaref'] + ':' + f1['aaalt']
    f1['mutation_ref'] = f1['gene']+';chr'+f1['chr'] + ';' + f1['start'] + ';'+f1['start'] +';'+ f1['ref'] + '>' + f1['alt']
    f1.to_csv(outpath + 'docm3.2_NS_index.tsv', sep='\t', index=False)


def clean(): 
    f = pd.read_csv(outpath + 'docm3.2_NS_index.tsv', sep='\t', dtype=str)
    need = f[base].drop_duplicates()
    need['source'] = 'DoCM'
    need.to_csv(outpath + "docm3.2_ns_clean.tsv",sep='\t',index=False)


def main():
    docmprocess1()
    docmprocess2()
    docmprocess3()
    docmprocess4()
    docmprocess5()
    clean()

if __name__ == '__main__':
	main()
