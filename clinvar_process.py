#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-07 11:05:22
Last Modified time: 2019-09-07 11:05:22
Description: Extract non-synonymous mutations in the clinvar dataset.
'''

import re
import pandas as pd


inpath = "/data/wangs/MutaAnalysis/clinvar/"
outpath = '/data/wangs/MutaAnalysis/clinvar/output/'


map = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
       "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
       "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
       "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
       "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
       "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
       "PTR": "X", "XLE": "X", "XAA": "X", "TER": "X"}


# Remove duplicate lines of the original file
def clinvar1():
    f1 = open(inpath + "variant_summary.txt", "r")
    f2 = open(outpath + "variant_summary_1.tsv", "w")
    a = set()
    for line in f1:
        if line not in a:
            f2.write(line)
            a.add(line)
    f2.close()


def clinvar2():  # output all variants of single nucleotide variant
    input_txt = open(outpath + "variant_summary_1.tsv", "r")
    output_txt = open(outpath + "filter_clinvar_1.tsv", "w")
    output_txt.write(input_txt.readline())
    for line in input_txt:
        temp = line.split('\t')  # 分割每一行
        if temp[1] != "single nucleotide variant":
            continue
        else:
            if temp[16] == "":  # temp[16] = "GRCh37" or "GRCh38"
                continue
            if temp[19] != temp[20]:  # strart position is different with stop position
            	continue
            if (temp[21] == "na") or (temp[22] == "na"):  # temp[21] = 'ReferenceAllele', temp[22] = 'AlternateAllele', although it's single nucleotide variant, but don't have allele information
            	continue
            if ("fs" in temp[2]) or ("del" in temp[2]) or ("ins" in temp[2]):
                continue
            else:
                output_txt.write(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
                    "\t%s"
                    "\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (temp[0], temp[1], temp[2], temp[3],
                                                      temp[4], temp[5], temp[6], temp[7], temp[8], temp[9],
                                                      temp[10], temp[11], temp[12], temp[13], temp[14],
                                                      temp[15],
                                                      temp[16], temp[17], temp[18], temp[19], temp[20],
                                                      temp[21],
                                                      temp[22], temp[23], temp[24], temp[25], temp[26],
                                                      temp[27], temp[28], temp[29]))
    input_txt.close()
    output_txt.close()


def clinvar3():  # delete "nsv/esv [dbVar]","Guidelines" column
    f1 = pd.read_table(outpath + 'filter_clinvar_1.tsv')
    f2 = f1.drop(['nsv/esv [dbVar]', 'Guidelines'], axis = 1)
    f2.to_csv(outpath + "filter_clinvar_2.tsv", sep = '\t', index = False)


def clinvar4():  # divide into GRCh37 and GRCh38
    f1 = pd.read_table(outpath + 'filter_clinvar_2.tsv')
    f2 = f1[f1['Assembly'] == 'GRCh37']
    f2.to_csv(outpath + 'ClinvarGRCh37.tsv', sep='\t', index=False)
    f3 = f1[f1['Assembly'] == 'GRCh38']
    f3.to_csv(outpath + 'ClinvarGRCh38.tsv', sep='\t', index=False)


def clinvar5():
    infile = open(outpath + 'ClinvarGRCh37.tsv', 'r')
    outfile1 = open(outpath + "ClinvarGRCh37_NS.tsv", "w")
    outfile2 = open(outpath + "ClinvarGRCh37_Silent.tsv", "w")
    for line in infile:
        temp = line.split('\t')
        if temp[0] == "#AlleleID":
            temp.insert(3, 'aaref')
            temp.insert(4, "aaalt")
            temp.insert(5, 'aachange')
            outfile1.write('\t'.join(temp))
            outfile2.write(line)
        if '?' in temp[2]:  # Tyr?His
            continue
        if '(p.=)' in temp[2]:  # NM_000044.4(AR):c.639G>A (p.=), don't have aa change
            continue
        else:
            a = temp[2].split(' ')
            if len(a) == 1:
                continue
            else:
                aa = a[1].strip(')').split('.')[1]
                pos = re.findall('\d+', aa)
                refaa = re.split('\d+', aa)[0]
                altaa = re.split('\d+', aa)[1]
                refaa = refaa.upper()
                altaa = altaa.upper()
                if (refaa == altaa) or ('=' in altaa):  # (p.Cys33Cys=), (p.Arg268=)
                    outfile2.write(line)
                else:
                    if len(refaa) == 1 and altaa == "*":  # p.Q1702*
                        temp.insert(3, refaa)
                        temp.insert(4, "X")
                        temp.insert(5, aa)
                    if len(refaa) == 1 and len(altaa) == 1 and altaa != "*":  # p.P672L
                        temp.insert(3, refaa)
                        temp.insert(4, altaa)
                        temp.insert(5, aa)
                    if len(refaa) == 3 and altaa == "*":
                        temp.insert(3, map[refaa])
                        temp.insert(4, "X")
                        temp.insert(5, '{}{}{}'.format(map[refaa], pos[0], altaa))
                    if len(refaa) == 3 and len(altaa) == 3:
                        temp.insert(3, map[refaa])
                        temp.insert(4, map[altaa])
                        if altaa == 'TER':
                            temp.insert(5, '{}{}{}'.format(map[refaa], pos[0], '*'))
                        else:
                            temp.insert(5, '{}{}{}'.format(map[refaa], pos[0], map[altaa]))
                    outfile1.write('\t'.join(temp))
    infile.close()
    outfile1.close()
    outfile2.close()


def clinvar6():
    f1 = pd.read_table(outpath + 'ClinvarGRCh37_NS.tsv', dtype=str)
    f1.insert(6,'aapos',f1['aachange'].str.findall('\d+').str[0])
    f1['mutationtype'] = ''
    f1.loc[f1['aaalt'] == 'X', 'mutationtype'] = 'nonsense'
    f1.loc[f1['aaalt'] != 'X', 'mutationtype'] = 'missense'
    f1['transcriptid'] = f1['Name'].str.split('(').str[0]
    f1['index'] = f1['Chromosome']+':'+f1['Start']+':'+f1['ReferenceAllele']+':'+f1['AlternateAllele']+':'+f1['aaref']+':'+f1['aaalt']
    f1.to_csv(outpath + 'ClinvarGRCh37_NS_index.tsv', sep='\t', index=False)


if __name__ == '__main__':
    clinvar1()
    clinvar2()
    clinvar3()
    clinvar4()
    clinvar5()
    clinvar6()
