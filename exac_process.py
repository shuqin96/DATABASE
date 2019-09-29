#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Time: 2019-08-27 10:06:08
Description: Processing ExAC non-TCGA version1.0 data.Select single nucleotide variants.
	We identified the variants leading to amino acid substitutions (AASs) by using the annotations from the Variant Effect Predictor (VEP) included in the downloaded VCF file.
'''

import os
import gc
import re
import pandas as pd
from itertools import cycle

path = "/data1/wangs/CanDriver/data/ExAC/data/"
outpath = "/data1/wangs/CanDriver/data/ExAC/output/"

##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP.
# Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position
# |Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|
# HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|
# EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|
# ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral

# CSQ=C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs554760071|1|880|-1||SNV|1
#|HGNC|38034||||||||||||||C:0.0020|C:0|C:0.0068|C:0.0014|C:0|C:0|C:0|||C:0|C:5.630e-05|C:0.003195|C:0.0001704|C:0|C:0|C:0|C:0|||||||||||||GGC|.


map = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
       "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
       "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
       "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
       "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
       "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
       "PTR": "X", "XLE": "X", "XAA": "X", "TER": "X", '*': 'X'}
base = ['gene','chr','start','ref','alt','aaref','aaalt','aapos','aachange','mutationtype','key','index','mutation_ref']


def filter_pass():
	os.chdir(('%s') % (path))
	os.system("java -jar /data/wangs/webserve/snpEff/SnpSift.jar filter '(FILTER = 'PASS')' ExAC_nonTCGA.r1.sites.vep.vcf > ExAC_nonTCGA.r1.sites.filter_PASS.vcf")


def excet_csq():
	os.chdir(('%s') % (path))
	os.system("java -jar /data/wangs/webserve/snpEff/SnpSift.jar extractFields -s ',' -e '.' ExAC_nonTCGA.r1.sites.filter_PASS.vcf CHROM POS REF ALT AF[*] CSQ > /data1/wangs/CanDriver/data/ExAC/output/ExAC_nonTCGA.r1.sites.filter_PASS_need_cols.tsv")


def select_single_nucleotide():
	f = pd.read_csv(outpath + "ExAC_nonTCGA.r1.sites.filter_PASS_need_cols.tsv",sep='\t',dtype=str)
	snp = f[f['REF'].isin(['A','T','G','C'])]
	snp.to_csv(outpath + "ExAC_nonTCGA.r1.sites.filter_PASS_snp.tsv",sep='\t',index=0)


def split_one_line():
    infile = open(outpath + "ExAC_nonTCGA.r1.sites.filter_PASS_snp.tsv", 'r')
    outfile = open(outpath + "ExAC_nonTCGA.r1.sites.filter_PASS_snp_split.tsv",'w')
    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('chr','start','ref','alt','af','AN_Adj','gene','mutationtype','aachange','aapos','exac_transcriptid'))
    infile.next()
    for line in infile:
        temp = line.strip('').split('\t')
        alt = temp[3]
        af = temp[4]
        csq = temp[6]
        if ',' not in alt: # only have one alternate allele
            csql = csq.split(',') # some variant have more than one effect
            for xinxi in csql:
                type = xinxi.split('|')[1]
                gene = xinxi.split('|')[3]
                enst = xinxi.split('|')[10]
                aapos = xinxi.split('|')[14]
                aa = xinxi.split('|')[15]
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],gene,type,aa,aapos,enst))
        else:  # There are two alts, the order of the amino acid is not necessarily the same as the previous one, so you have to match the corresponding ALT
            altlist = alt.split(',')
            aflist = af.split(',')
            csql = csq.split(',')
            for altalt in altlist:
                id = altlist.index(altalt) # find the position of ALT
                for xinxi in csql:
                    allele = xinxi.split('|')[0]
                    type = xinxi.split('|')[1]
                    gene = xinxi.split('|')[3]
                    enst = xinxi.split('|')[10]
                    aapos = xinxi.split('|')[14]
                    aa = xinxi.split('|')[15]
                    if altalt == allele: # Make sure that the corresponding ALT is found
                        outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (temp[0],temp[1],temp[2],altalt,aflist[id],temp[5],gene,type,aa,aapos,enst))
                    else:
                        continue
    outfile.close()


def coding_snp_muta():
	f = pd.read_table(outpath + "ExAC_nonTCGA.r1.sites.filter_PASS_snp_split.tsv",dtype=str)
	f = f.drop_duplicates()  
	base = ['A','T','G','C']
	snp = f[f['ref'].isin(base) & (f['alt'].isin(base))] 
	# missense_variant&splice_region_variant 
	cod = snp[snp['mutationtype'].str.contains('missense_variant|stop_gained')] 
	cod['ENST'] = cod['exac_transcriptid'].str.split(':').str[0]
	cod['cDNA'] = cod['exac_transcriptid'].str.split(':').str[1]
	cod['aaref'] = cod['aachange'].str.split('/').str[0]
	cod['aaalt'] = cod['aachange'].str.split('/').str[1]
	cod['aachange'] = cod['aaref']+cod['aapos']+cod['aaalt']
	cod['aaalt'].replace('*','X',inplace=True)
	cod['mutationtype'].replace({'missense_variant':'missense','stop_gained':'nonsense'},inplace=True)
	# cod.loc[cod['mutationtype'].str.contains('missense_variant'), 'mutationtype'] = 'missense'
	# cod.loc[cod['mutationtype'].str.contains('stop_gained'), 'mutationtype'] = 'nonsense'
	cod.drop(['exac_transcriptid'],axis=1,inplace=True)
	cod.to_csv(outpath + "ExAC_nonTCGA.r1.sites.ns.tsv",sep='\t',index=False)


def add_index():
    f = pd.read_table(outpath + "ExAC_nonTCGA.r1.sites.ns.tsv", dtype=str) # 6310383
    f['mutation_ref'] = f['gene']+';chr'+f['chr']+';'+f['start']+';'+f['start']+';'+f['ref']+'>'+f['alt']
    f['index'] = f['chr']+':'+f['start']+':'+f['ref']+':'+f['alt']+':'+f['aaref']+':'+f['aaalt']
    f['key'] = f['gene']+':'+f['aachange']
    f.to_csv(outpath + "ExAC_nonTCGA.r1.sites.ns.index.tsv",sep='\t',index=False)


def clean():
	f = pd.read_csv(outpath + "ExAC_nonTCGA.r1.sites.ns.index.tsv",sep='\t',dtype=str)
	need = f[base+['af']]
	need = need.drop_duplicates()
	need.to_csv(outpath + "ExAC_nonTCGA.r1.sites.ns.index.clean.tsv",sep='\t',index=0)


## 计算各种分数,这里计算时用的是最开始生成的文件，没有判断是否为编码，所以比较多
def get_score():
	os.chdir('/data/wangs/MutaAnalysis/')
	os.system("python data_process.py -p %s -d ExAC_nonTCGA -i ExAC_nonTCGA.r1.sites.ns.index.clean.tsv" % (outpath))


## candra数据量多时计算地太慢，所以，我直接用之前计算好的，把之前没有的部分重新计算即可
def candra_new():
	f1 = pd.read_csv(outpath + "ExAC_nonTCGA_candra_input.txt",sep='\t',dtype=str,names=['chr','start','ref','alt','chain'])
	f2 = pd.read_csv("/data1/wangs/CanDriver/data/ExAC/output_early/ExAC_nonTCGA_candra_input.txt",sep='\t',dtype=str)
	f1['index']=f1['chr']+':'+f1['start']+':'+f1['ref']+':'+f1['alt']
	f2['index']=f2['chr']+':'+f2['start']+':'+f2['ref']+':'+f2['alt']
	f1index = set(f1['index'])
	f2index = set(f2['index'])
	new=f1index-f2index
	out = f1[f1['index'].isin(list(new))]
	out.drop(['index'],axis=1,inplace=True)
	out.to_csv(outpath + "ExAC_nonTCGA_candra_duo_input.txt",sep='\t',index=0)
	# get score
	os.chdir('%s' % (outpath))
	os.system("perl /data/wangs/webserve/CanDrA.v+/open_candra.pl GENERAL ./ExAC_nonTCGA_candra_duo_input.txt > ExAC_nonTCGA_candra_duo_output.tsv")


def candra_out_combine():
	f1 = pd.read_csv(outpath + "ExAC_nonTCGA_candra_duo_output.tsv",sep='\t',dtype=str)
	f2 = pd.read_csv("/data1/wangs/CanDriver/data/ExAC/output_early/ExAC_nonTCGA_candra_output.tsv",sep='\t',dtype=str)
	fw = pd.concat([f1,f2],axis=0)
	fw.to_csv(outpath +  "ExAC_nonTCGA_candra_output.tsv",sep='\t',index=0)


if __name__ == '__main__':
	# select_single_nucleotide()
	# split_one_line()
	# coding_snp_muta()
	# add_index()
	# clean()
	# get_score()
	# candra_new()
	# candra_out_combine()
	get_score() # 从merge_all_score开始
