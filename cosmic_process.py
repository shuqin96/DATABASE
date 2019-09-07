#!/usr/bin/python
# -*- coding: UTF-8 -*-

'''
@Author: Wangshuqin
@Date:   2019-03-18 18:55:27
@Last Modified by:   WangShuqin
@Last Modified time: 2019-06-26 09:30:57
@Description: Extract non-synonymous mutations in the COSMIC v87 dataset. Add CancerGeneCensus informations and add ICGC cancer type by NCI code. 
'''

import pandas as pd
import gc,re,os

inpath = "/data1/wangs/CanDriver/data/COSMIC/v87/"
outpath = "/data1/wangs/CanDriver/data/COSMIC/v87/output/"

# For cDNA to gDNA correspondence
map = {'A': 'T',
		'T': 'A',
		'G': 'C',
		'C': 'G'}

# Select single nucleotide variant
def cosmic_single_nucleotide():
	f1 = pd.read_csv(inpath + 'CosmicGenomeScreensMutantExport.tsv',sep='\t',dtype=str)
	f1['variant'] = f1['Mutation AA'].str.split('.').str[1]  # p.R184C
	f1['Chr'] = f1['Mutation genome position'].str.split(':').str[0]  # 12:25398285-25398285
	# 人类有23对染色体，22对常染色体和1对性染色体，cosmic中23号染色体表示的是X染色体，24表示的是Y, 但是其中还有25号染色体，它的基因名中都有MT，所以我觉得是表示线粒体DNA
	f1['Chr'] = f1['Chr'].replace('23', 'X')
	f1['Chr'] = f1['Chr'].replace('24', 'Y')
	f1['Start'] = f1['Mutation genome position'].str.split('[:,-]').str.get(1)
	f1['Stop'] = f1['Mutation genome position'].str.split('[:,-]').str.get(2)
	f1['cDNA_ref'] = f1['Mutation CDS'].str.get(-3)  # temp[17] = "c.711C>T"
	f1['cDNA_alt'] = f1['Mutation CDS'].str.get(-1)
	f1['Ref_AA'] = f1['variant'].str.split('\d+').str[0]
	f1['Alt_AA'] = f1['variant'].str.split('\d+').str[1]
	f1['Ref_AA'] = f1['Ref_AA'].replace('*','X')
	f1['Alt_AA'] = f1['Alt_AA'].replace('*','X')
	f1['aapos'] = f1['variant'].str.findall('\d+').str[0]
	# Remove mutations no genomic location information, no alternate amino acid information and on chromosome 25, Nonstop extension（*>*）这种不要
	f2 = f1[(f1['Mutation genome position'].notnull()) & f1['Mutation Description'].str.contains('Substitution') & (f1['Start'] == f1['Stop']) & (f1['Chr'] != '25')]
	# Variant Type, Substitution - Missense变为Missense
	f2['Mutation Description'] = f2['Mutation Description'].str.split('-').str[1].str.strip() 
	# Delete unuseful column
	f2.drop(['Pubmed_PMID'], axis=1, inplace=True)
	# 由于删除了一列，所以剩下的会有重复的
	fw = f2.drop_duplicates()
	fw.to_csv(outpath + 'CosmicSingleNucleotideMutant.tsv', sep='\t', index=False)


def add_ref_alt():
	# 在整个cosmic 负链上的ref,alt allele转变过来
	f = pd.read_csv(outpath + "CosmicSingleNucleotideMutant.tsv", dtype=str,sep='\t')
	f.loc[f['Mutation strand'] == '-', 'ref'] = f.loc[f['Mutation strand'] == '-', 'cDNA_ref'].map(map)
	f.loc[f['Mutation strand'] == '-', 'alt'] = f.loc[f['Mutation strand'] == '-', 'cDNA_alt'].map(map)
	f.loc[f['Mutation strand'] == '+', 'ref'] = f.loc[f['Mutation strand'] == '+', 'cDNA_ref']
	f.loc[f['Mutation strand'] == '+', 'alt'] = f.loc[f['Mutation strand'] == '+', 'cDNA_alt']
	f.to_csv(outpath + "CosmicSingleNucleotideMutantAddallele.tsv",sep='\t',index=False)


def drop_gene():
	# The gene name have many : ENSG00000257591\Q96QE0_HUMAN\NRG1_ENST00000521670\FBXW7_NM_018315_2\TERT_NM_198255.1\NM_198455_2\NP_001073948_1\NR_003148_2
	f = pd.read_csv(outpath + "CosmicSingleNucleotideMutantAddallele.tsv", dtype=str, sep='\t')
	fg = set(f['Gene name'])
	nm = f[f['Gene name'].str.startswith('NM_')|f['Gene name'].str.startswith('NP_')|f['Gene name'].str.startswith('NR_')]
	nm['Gene'] = nm['Gene name']
	nmg = set(nm['Gene'])
	hum = f[f['Gene name'].str.contains('HUMAN')]
	hum['Gene'] = hum['Gene name']
	humg = set(hum['Gene'])
	og = fg-nmg-humg
	other = f[f['Gene name'].isin(list(og))]
	other['Gene'] = other['Gene name'].str.split('_').str[0]
	fw = pd.concat([nm,hum,other],axis=0)
	fw.to_csv(outpath + "CosmicSingleNucleotideMutantGene.tsv",sep='\t',index=False)
	# we only select mutations whose gene name is normal
	other.to_csv(outpath + "CosmicSingleNucleotideMutantGeneName.tsv",sep='\t',index=False)


# cosmic数据集中，gene name中有ENST, 另外有一列也是transcript编号，会有不一致的情况，这种的就将其去除
def drop_error():
	cosmic = pd.read_csv(outpath + "CosmicSingleNucleotideMutantGeneName.tsv", dtype=str, sep='\t')
	# 因为这一列里既有NM，又有ENST，所以分开处理
	cosmic.loc[cosmic['Accession Number'].str.startswith('ENST'), 'accession_number'] = cosmic.loc[cosmic['Accession Number'].str.startswith('ENST'),'Accession Number'].str.split('_').str[0]
	cosmic.loc[~cosmic['Accession Number'].str.startswith('ENST'), 'accession_number'] = cosmic.loc[~cosmic['Accession Number'].str.startswith('ENST'), 'Accession Number']
	# gene名里有enst的，gene_enst就是gene里的;gene名里没enst的gene_enst就是Assession number里提取出的最终的
	cosmic.loc[cosmic['Gene name'].str.contains("ENST"), 'gene_enst'] = cosmic.loc[cosmic['Gene name'].str.contains("ENST"), 'Gene name'].str.split('_').str[1]
	cosmic.loc[~cosmic['Gene name'].str.contains("ENST"), 'gene_enst'] = cosmic.loc[~cosmic['Gene name'].str.contains("ENST"), 'accession_number']
	exact = cosmic[cosmic['gene_enst'] == cosmic['accession_number']]
	exact.drop(['Gene name','Accession Number','gene_enst'], axis=1, inplace=True)
	exact.to_csv(outpath + "CosmicSingleNucleotideMutantClean.tsv",sep='\t',index=False)
	

# Add informations from cancer census gene, such as ‘TSG/Oncogene’...
def add_gene_census_info():
	f = pd.read_csv(outpath + "CosmicSingleNucleotideMutantClean.tsv", dtype=str, sep='\t')
	census = pd.read_csv(inpath + "cancer_gene_census.csv", dtype=str)
	fw = f.merge(census[['Gene Symbol', 'Tier','Molecular Genetics', 'Role in Cancer']],
                   how='left',left_on=['Gene'],right_on=['Gene Symbol'])
	fw.drop(['Gene Symbol'],axis=1,inplace=True)
	fw.to_csv(outpath + "CosmicSNVaddCensusINFO.tsv",sep='\t',index=False)


def somatic():
	f = pd.read_csv(outpath + "CosmicSNVaddCensusINFO.tsv", dtype=str, sep='\t')
	somatic = f[f['Mutation somatic status'] == 'Confirmed somatic variant']
	# Add index
	somatic['mutation_ref'] = somatic['Gene']+';chr'+somatic['Chr']+';'+somatic['Start']+';'+somatic['Stop']+';'+somatic['ref']+'>'+somatic['alt']
	somatic['index'] = somatic['Chr']+':'+somatic['Start']+':'+somatic['ref']+':'+somatic['alt']+':'+somatic['Ref_AA']+':'+somatic['Alt_AA']
	somatic.to_csv(outpath + 'CosmicSomaticMutant.tsv', sep='\t', index=False)


def add_cancer_type():
	f1 = pd.read_csv(outpath + "CosmicSomaticMutant.tsv", dtype=str, sep='\t')
	f1['ID_STUDY'] = f1['ID_STUDY'].str.split('.').str[0]
	# 1.add nci code
	sample = pd.read_csv(inpath + "CosmicSample.tsv", sep='\t', dtype=str)
	fw1 = f1.merge(sample[['sample id', 'sample_name', 'NCI code']], how='left', left_on=['ID_sample','Sample name'],right_on=['sample id', 'sample_name'])
	del sample,f1
	# 2.cosmic网页上有关cosmic同icgc cancer type的匹配
	cancer = pd.read_excel("/data1/suntt/CanDriver/Data/COSMIC_study_id.xlsx",dtype=str)
	fw2 = fw1.merge(cancer[['Study','Tumour Type','Study Id']], how='left', left_on=['ID_STUDY'], right_on=['Study Id'])
	fw3 = fw2[fw2['Study'].notnull()] # 能匹配到ICGC的部分
	non = fw2[fw2['Study'].isnull()] # 不能匹配到ICGC的部分
	del fw1,fw2
	# 没有与icgc匹配到的部分，用nci code对应的cancer名称
	nci = pd.read_csv(inpath + "cosmic_nci_code_cancer.txt", sep='\t',dtype=str)
	fw4 = non.merge(nci, how='left', left_on=['NCI code'], right_on=['nci_code'])
	del nci
	gc.collect()
	# combine
	fw = pd.concat([fw3,fw4], axis=0)
	fw.drop(['nci_code','Study Id','sample id','sample_name'],axis=1,inplace=True)
	fw.rename(columns={'Study':'icgc_project_code'},inplace=True)
	fw.to_csv(outpath + "CosmicSomaticMutantCancerType.tsv",sep='\t',index=False)


## Filter the required data based on criteria
def filter_data():
	f1 = pd.read_csv(outpath + "CosmicSomaticMutantCancerType.tsv", dtype=str, sep='\t')
	# 1.samples which came from cell-lines, xenografts, or organoid cultures were excluded.
	f2 = f1[~f1['Sample Type'].isin(['cell-line', 'xenograft', 'organoid culture'])]
	# 2.mutations which were flagged as SNPs were excluded.
	f3 = f2[f2['SNP'] != 'y']  # Mutations no this information is also retained, otherwise there will be too much mutations are removed
	f3.to_csv(outpath + "CosmicSomaticMutantFilter.tsv", sep='\t', index=False) 
	# 3.For analysis comparing oncogenes and tumor suppressor genes (TSG), genes classified as only fusion genes or those with both oncogenic and TSG activities were not used.
	f4 = f3[(~f3['Role in Cancer'].isin(['fusion', 'oncogene, TSG', 'oncogene, TSG, fusion'])) & (f3['Role in Cancer'].notnull())]
	f4.to_csv(outpath + "CosmicSomaticMutantFilterGene.tsv", sep='\t', index=False)
	# print(f1.shape[0], f2.shape[0],f3.shape[0],f4.shape[0])


def main():
	cosmic_single_nucleotide()
	add_ref_alt()
	drop_error()
	add_gene_census_info()
	somatic()  # and add index
	add_cancer_type()
	filter_data()


if __name__ == '__main__':
	main()
