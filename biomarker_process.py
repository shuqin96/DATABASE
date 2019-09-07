#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-07 12:41:42
Last Modified time: 2019-09-07 12:41:42
Description: Extract non-synonymous mutations in the CGI Biomarker dataset.
'''

import pandas as pd

inpath = '/data/wangs/MutaAnalysis/Hpm/biomarker/'
outpath = '/data/wangs/MutaAnalysis/Hpm/biomarker/output/'

data = 'cgi_biomarkers'
base = ['gene','chr','start','ref','alt','aaref','aaalt','aapos','aachange','mutationtype','key','index','mutation_ref']

def snv():
	f1 = pd.read_table(inpath + 'cgi_biomarkers_per_variant.tsv')
	# 有氨基酸变异的信息
	f2 = f1[f1['individual_mutation'].notnull() & f1['gDNA'].notnull()]
	# single nucleotide variants
	# 虽然有些变体有氨基酸酸变异的信息，但他们可能是由delete or insert变化引起的, 所以得排除
	f3 = f2[~(f2['gDNA'].str.contains('del'))]
	f3.rename(columns={'Gene':'gene'},inplace=True)
	f3.to_csv(outpath + 'cgi_biomarkers_snv.tsv', sep='\t', index=False)


def add_info():
	f = pd.read_table(outpath + 'cgi_biomarkers_snv.tsv')
	f['mutationtype'] = f['info'].str.split('[=;]').str[1]
	f['chr'] = f['gDNA'].str.split(':').str[0].str.strip('chr')
	f['start'] = f['gDNA'].str.findall('\d{3,}').str[0]
	f['ref'] = f['gDNA'].str[-3]
	f['alt'] = f['gDNA'].str[-1]
	f['aaref'], f['aaalt'] = f['individual_mutation'].str.split(':').str[1].str.split('\d+').str
	f['aaalt'].replace('*', 'X', inplace=True)
	f['aachange'] = f['individual_mutation'].str.split(':').str[1]
	f['aapos'] = f['aachange'].str.findall('\d+').str[0]
	ns = f[f['mutationtype'].isin(['Missense','Nonsense'])]
	ns['mutationtype'] = ns['mutationtype'].str.lower()
	ns.to_csv(outpath + 'cgi_biomarkers_ns_add_info.tsv', sep='\t', index=False)


# cgi_biomarkers.tsv in which variants that constitute biomarkers of the same type of response response to the same drug in a given cancer type(s) are grouped in the same row.
# cgi_biomarkers_per_variant.tsv each variant is contained in a separate row and --for point mutations-- the genomic coordinates are included.
def group_by_drug():
	f1 = pd.read_table(outpath + 'cgi_biomarkers_ns_add_info.tsv')
	group = f1.groupby('individual_mutation')
	col = ['Association', 'Drug', 'Drug family', 'Drug full name', 'Drug status', 'Evidence level']
	group1 = group[col].agg(lambda x:'|'.join(x)).reset_index()
	fw = f1.drop(col, axis=1).merge(group1, how='left', on=['individual_mutation'])
	uni = fw.drop_duplicates('individual_mutation')
	uni.to_csv(outpath + 'cgi_biomarkers_ns_unique.tsv', sep='\t', index=False)


def add_index():
	f = pd.read_table(outpath + 'cgi_biomarkers_ns_unique.tsv', dtype=str)
	f['key'] = f['gene']+':'+f['aachange']
	f['index'] = f['chr'] + ':' + f['start'] + ':' + f['ref'] + ':' + f['alt'] + ':' + f['aaref'] + ':' + f['aaalt']
	f['mutation_ref'] = f['gene']+';chr'+f['chr'] + ';' + f['start'] + ';'+f['start'] +';'+ f['ref'] + '>' + f['alt']
	f.to_csv(outpath + 'cgi_biomarkers_ns_index.tsv', sep='\t', index=False)


def clean():
	f=pd.read_csv(outpath + 'cgi_biomarkers_ns_index.tsv', sep='\t', dtype=str)
	need = f[base+['Association']].drop_duplicates()
	need['source'] = 'Biomarker'
	need.to_csv(outpath + data + "_ns_clean.tsv",sep='\t',index=False)

def damage():
    f = pd.read_csv(outpath +data +"_ns_clean.tsv",sep='\t',dtype=str)
    f = f[~f['Association'].str.contains('No Responsive')]
    f.to_csv(outpath+data+"_ns_clean_damage_mutation.tsv",sep='\t',index=False)
      

def main():
	snv()
	add_info()
	group_by_drug()
	add_index()
	clean()
	damage()

if __name__ == '__main__':
	main()
