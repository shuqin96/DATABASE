#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-03 20:23:32
Last Modified time: 2019-09-03 20:23:32
Description: A mutation of a genome position can have multiple scores.
			SIFT_score、Polyphen2_HDIV_score、Polyphen2_HVAR_score、MutationAssessor_score、PROVEAN_score、VEST4_score have multiple scores, corresponding to Ensembl_proteinid/Ensembl_transcriptid/Uniprot_acc/Uniprot_entry.
			FATHMM_score have multiple scores corresponding to Ensembl_proteinid, but there are special condition that need to be handled separately.
			MutationTaster_score have multiple scores by its MutationTaster_AAE(amino acid change).
			There are Amino acid change used for MutPred_score calculation.
			The script is written according to the hg38 location. If your file is hg19, you can replace the 'pos(1-based)' column with the 'hg19_pos(1-based)' column.
'''


import pandas as pd
import sys, getopt

data = ''
outpath = ''
# set up input name and outpath for each job
opts, args = getopt.getopt(sys.argv[1:], "d:p:", ['data=', 'path='])
for name, value in opts:
    if name in ('-d', '--data'): # 命令行参数-d, --data指定的是该脚本中data的值
        data = value
    if name in ('-p', '--path'): # 命令行参数-p, --path指定的是该脚本中outpath的值
        outpath = value

data = data
outpath = outpath


# Select the required columns according to your needs
def dbnsfp_clean():
    dbnsfp = pd.read_csv(outpath + data + '_dbnsfp_output.tsv', sep='\t', dtype=str)
	dbnsfp = dbnsfp[['#chr','pos(1-based)','ref','alt','aaref','aaalt','aapos']]
	dbnsfp.to_csv(outpath + data + '_dbnsfp_output_clean.tsv',sep='\t',index=0)


# 以下方法的score是按aapos(ENST,ENSP,UniprotID)对应分开的,但是有个例外，FATHMM score虽然说是按ENSP对应的，但是当它多个对应的都是'.'时，它只展示一个'.',这样就不能用以下程序对应，得单独处理
def dbnsfp1():
	df = pd.read_csv(outpath + data + '_dbnsfp_output_clean.tsv',sep='\t',dtype=str)
	# These columns need to be processed separately, so remove them first.
	df = df.drop(['MutationTaster_score','MutationTaster_pred','MutationTaster_AAE','FATHMM_score','FATHMM_pred','FATHMM_converted_rankscore','MutPred_score','MutPred_rankscore','MutPred_protID','MutPred_AAchange','MutPred_Top5features'],axis=1)
	aapos = df['aapos'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('aapos')
	unp = df['Uniprot_acc'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Uniprot_acc')
	sift = df['SIFT_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('SIFT_score')
	siftp = df['SIFT_pred'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('SIFT_pred')
	pph2 = df['Polyphen2_HDIV_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Polyphen2_HDIV_score')
	pph2p = df['Polyphen2_HDIV_pred'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Polyphen2_HDIV_pred')
	pph2a = df['Polyphen2_HVAR_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Polyphen2_HVAR_score')
	pph2ap = df['Polyphen2_HVAR_pred'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Polyphen2_HVAR_pred')
	ma = df['MutationAssessor_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('MutationAssessor_score')
	map = df['MutationAssessor_pred'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('MutationAssessor_pred')
	prv = df['PROVEAN_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('PROVEAN_score')
	prvp = df['PROVEAN_pred'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('PROVEAN_pred')
	vest = df['VEST4_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('VEST4_score')
	fw = pd.concat([aapos,unp,sift,siftp,pph2,pph2p,pph2a,pph2ap,ma,map,prv,prvp,vest],axis=1)
	# Merged with the original file, so the output is the entire file, not just the columns.
	dfw=df.drop(['aapos','Uniprot_acc','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred',
		'MutationAssessor_score','MutationAssessor_pred','PROVEAN_score','PROVEAN_pred','VEST4_score'], axis=1).join(fw)
	# After this separation, there are multiple rows of the same mutation, so remove the duplicate and keep the first row.
	uni = dfw.drop_duplicates()
	uni.to_csv(outpath + data + "_dbnsfp_output_split.tsv",sep='\t',index=0)


def dbnsfp_mutationtaster():
	df = pd.read_csv(outpath + data + '_dbnsfp_output.tsv',sep='\t',dtype=str)
	# select MutationTaster columns
	df = df[['#chr','pos(1-based)','ref','alt','aaref','aaalt','MutationTaster_score','MutationTaster_pred','MutationTaster_AAE']]
	mt = df['MutationTaster_score'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('MutationTaster_score')
	mtp = df['MutationTaster_pred'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('MutationTaster_pred')
	mta = df['MutationTaster_AAE'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('MutationTaster_AAE')
	fw = pd.concat([mt,mtp,mta],axis=1)
	dfw=df.drop(['MutationTaster_score','MutationTaster_pred','MutationTaster_AAE'], axis=1).join(fw)
	uni = dfw.drop_duplicates()
	uni.to_csv(outpath + data + "_dbnsfp_output_mutationtaster.tsv",sep='\t',index=0)


def dbnsfp_fathmm():
	df = pd.read_csv(outpath + data + '_dbnsfp_output.tsv',sep='\t',dtype=str)
	fw = open(outpath + data + "_dbnsfp_output_fathmm.tsv",'w')
	fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('#chr','pos(1-based)','ref','alt','aaref','aaalt','aapos','FATHMM_score','FATHMM_pred','FATHMM_converted_rankscore'))
	# select FATHMM columns
	df = df[['#chr','pos(1-based)','ref','alt','aaref','aaalt','aapos','FATHMM_score','FATHMM_pred','FATHMM_converted_rankscore']]
	for index, row in df.iterrows():
		aapos_list = row['aapos'].split(';')
		famlist = row['FATHMM_score'].split(';')
		famplist = row['FATHMM_pred'].split(';')
		for aapos,fam,famp in zip(aapos_list,cycle(famlist),cycle(famplist)):
			fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (row['#chr'],row['pos(1-based)'],row['ref'],row['alt'],row['aaref'],row['aaalt'],aapos,fam,famp,row['FATHMM_converted_rankscore']))


def dbnsfp_mutpred(): 
	df = pd.read_csv(outpath + data + '_dbnsfp_output.tsv',sep='\t',dtype=str)
	# select MutPred columns
	df = df[['#chr','pos(1-based)','ref','alt','aaref','aaalt','MutPred_score','MutPred_rankscore','MutPred_protID','MutPred_AAchange']]
	df = df.drop_duplicates()
	df.to_csv(outpath + data + "_dbnsfp_output_mutpred.tsv",sep='\t',index=0)


# infile:This is the file you need for dbNSFP scores. 
def merge():
	f = pd.read_csv(outpath + infile,sep='\t',dtype=str)
	# merge methods correspond to aapos
	unp = pd.read_csv(outpath + data + "_dbnsfp_output_split.tsv",sep='\t',dtype=str) 
	fw = f.merge(unp, how='left', left_on=[chrcol, poscol, refcol, altcol, aarefcol, aaaltcol, aaposcol], right_on=['#chr','pos(1-based)','ref','alt','aaref','aaalt','aapos'])
	# merge FATHMM
	fathmm = pd.read_csv(outpath + data + "_dbnsfp_output_fathmm.tsv",sep='\t',dtype=str)
	uni = fathmm.drop_duplicates()
	fw = fw.merge(uni, how='left', left_on=[chrcol, poscol, refcol, altcol, aarefcol, aaaltcol, aaposcol], right_on=['#chr','pos(1-based)','ref','alt','aaref','aaalt','aapos'])
	# merge MutationTaster
	mt = pd.read_csv(outpath + data + "_dbnsfp_output_mutationtaster.tsv",sep='\t',dtype=str)
	fw = fw.merge(mt, how='left', left_on=[chrcol, poscol, refcol, altcol, aachangecol], right_on=['#chr','pos(1-based)','ref','alt','MutationTaster_AAE'])
	# merge mutpred
	mutpred = pd.read_csv(outpath + data + "_dbnsfp_output_mutpred.tsv",sep='\t',dtype=str)
	fw = fw.merge(mt, how='left', left_on=[chrcol, poscol, refcol, altcol, aachangecol], right_on=['#chr','pos(1-based)','ref','alt','MutPred_AAchange'])
	# Remove extra columns
	fw.drop(['#chr_x','pos(1-based)_x','ref_x','alt_x','MutationTaster_AAE','#chr_y','pos(1-based)_y','ref_y','alt_y','aaref_y','aaalt_y','aapos_y'],axis=1,inplace=True)
	fw.to_csv(outpath + data + "_merge_dbsnfp.tsv",sep='\t',index=0)


if __name__ == '__main__':
	dbnsfp1()
	dbnsfp_mutationtaster()
	dbnsfp_fathmm()
	dbnsfp_mutpred()
	merge()
