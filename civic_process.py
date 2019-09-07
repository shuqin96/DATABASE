#!/usr/bin/python
#-*- coding:utf-8 -*- 
'''
Author: Wangshuqin
Data: 2019-09-07 10:58:17
Last Modified time: 2019-09-07 10:58:17
Description: processing CIViC DATA
'''

import re
import pandas as pd
from itertools import cycle

inpath = "/data/wangs/MutaAnalysis/civic/"
outpath = '/data/wangs/MutaAnalysis/civic/output/'


def civic1():
    # output all information of single nucleotide variant
    input_txt = open(inpath + "nightly-ClinicalEvidenceSummaries.tsv", "r")
    output_txt = open(outpath + "civic_evidence_1.tsv", "w")
    output_txt.write(input_txt.readline())
    for line in input_txt:
        temp = line.split('\t')
        if ('FS' in temp[2]) or ('del' in temp[2]) or ('DEL' in temp[2]):
            continue
        if temp[19] == "":  # 排除没有位置信息的个体
            continue
        if temp[19] != temp[20]:  # 排除不是点突变的个体
            continue
        if temp[21] == '':
            continue
        elif re.search('[A-Z]{1}\d+[A-Z]{1}|[A-Z]\d+\*|\*\d+[A-Z]', temp[2]):
            # print temp[2]
            output_txt.write(line)
    input_txt.close()
    output_txt.close()


def civic2():  # delete unuseful columns
    f1 = pd.read_table(outpath + "civic_evidence_1.tsv")
    unuseful = ['pubmed_id', 'citation', 'chromosome2', 'start2', 'stop2', 'representative_transcript2', 'ensembl_version', 'evidence_civic_url', 'variant_civic_url', 'gene_civic_url']
    f2 = f1.drop(unuseful, axis=1)
    f2.to_csv(outpath + "civic_evidence_2.tsv", sep='\t', index=False)


# divide civic data into NS and Silent
def civic3():
    infile = open(outpath + "civic_evidence_2.tsv", "r")
    outfile1 = open(outpath + "CivicSilent.tsv", "w")
    outfile2 = open(outpath + "CivicNs.tsv", "w")
    outfile2.write(infile.readline())
    for line in infile:
        temp = line.split('\t')
        a = temp[2].split(' ')
        if "(" in temp[2]:  # A149T (c.445G>A)
            b = a[0]
            c = re.split('\d+', b)
            if c[0] == c[1]:
                outfile1.write(line)
            if c[0] != c[1]:
                outfile2.write(line)
        else:
            d = a[-1]
            e = re.split('\d+', d)
            if e[-1] == "":
                continue
            if len(e) == 1:  # TKD MUTATION
                continue
            else:
                if e[0] == e[1]:
                    outfile1.write(line)
                if e[0] != e[1]:
                    outfile2.write(line)
    outfile1.close()
    outfile2.close()


# civicNS insert refaa, altaa, aa position and variant columns
def civic4():
    infile = open(outpath + "CivicNs.tsv", "r")
    outfile = open(outpath + "CivicNs_add_aa.tsv", "w")
    for line in infile:
        temp = line.split('\t')
        if temp[2] == "variant":
            temp.insert(3, "ref_aa")
            temp.insert(4, "alt_aa")
            temp.insert(5,'aapos')
            temp.insert(6, "aachange")
        else:
            a = re.findall('[A-Z]\d+[A-Z]|\*\d+[A-Z]|[A-Z]\d+\*', temp[2])
            b = a[-1]
            c = re.split('\d+', b)
            f = re.findall('\d+', b)
            if c[0] == "*":
                temp.insert(3, 'X')
                temp.insert(4, c[1])
                temp.insert(5, f[0])
                temp.insert(6, b)
            if c[1] == '*':
                temp.insert(3, c[0])
                temp.insert(4, 'X')
                temp.insert(5, f[0])
                temp.insert(6, b)
            if '*' not in b:
                temp.insert(3, c[0])
                temp.insert(4, c[1])
                temp.insert(5, f[0])
                temp.insert(6, b)
        outfile.write('\t'.join(temp))


# Civic_add_aa insert another column mutation type
def civic5():
    infile = open(outpath + "CivicNs_add_aa.tsv", "r")
    outfile = open(outpath + "CivicNs_add_mutype.tsv", "w")
    for line in infile:
        temp = line.split('\t')
        if temp[4] == "alt_aa":
            temp.insert(7, "mutationtype")
        else:
            if temp[4] == "X":
                temp.insert(7, "nonsense")
            else:
                temp.insert(7, "missense")
        outfile.write('\t'.join(temp))
    outfile.close()


# Civic insert another collumn of pathogenic prediction,such as 'high','medium','low'
def civic6():
    infile = open(outpath + "CivicNs_add_mutype.tsv", "r")
    outfile = open(outpath + "CivicNs_add_prediction.tsv", "w")
    for line in infile:
        temp = line.split('\t')
        if temp[12] == "evidence_direction":
            temp.insert(15, "pathogenic_prediction")
        if temp[12] == "Does Not Support" or temp[12] == "":
            temp.insert(15, "nan")
        else:
            if temp[14] == "Pathogenic" or temp[14] == "Sensitivity" or temp[14] == "Poor Outcome":
                temp.insert(15, "High")
            if temp[14] == "Better Outcome" or temp[14] == "Positive":
                temp.insert(15, "Medium")
            if temp[14] == "Resistance or Non-Response" or temp[14] == "N/A" or temp[14] == "Uncertain Significance" \
                    or temp[14] == "Negative":
                temp.insert(15, "Low")
            if temp[14] == "":
                # print temp[12]
                temp.insert(15, "nan")
            if temp[14] == "Adverse Response":
                temp.insert(15, "nan")
        outfile.write('\t'.join(temp))
    outfile.close()


# Civic insert accession number column
def civic7():
    infile = open(outpath + "CivicNs_add_prediction.tsv", "r")
    outfile = open(outpath + "CivicNs_add_accessnum.tsv", "w")
    for line in infile:
        temp = line.split('\t')
        if temp[0] == 'gene':
            temp.insert(28, 'ensembl_number')
        else:
            a = temp[27].split('.')
            if a[0] != '':
                temp.insert(28, a[0])
            else:
                temp.insert(28, 'nan')
        outfile.write('\t'.join(temp))


# Civic insert index, '9:5073770:G:T:V:F'

def a(row):  # output chr, ref, alt, refaa, altaa
    return '{}:{}:{}:{}:{}:{}'.format(row[22], row[23], row[25], row[26], row[3], row[4])


def civic8():
    f1 = pd.read_table(outpath + "CivicNs_add_accessnum.tsv")
    f1['index'] = f1.apply(a, axis=1)
    f2 = f1.fillna('nan')
    f2.to_csv(outpath + "CivicNs_add_index.tsv", sep='\t', index=False)


def civic9():
    infile = open(outpath + "CivicNs_add_index.tsv", "r")
    outfile = open(outpath + "CivicNs_column6.tsv", "w")
    a = set()
    infile.next()
    for line in infile:
        temp = line.strip().split('\t')
        a.add(temp[-1])
    for i in a:
        b = i.split(':')
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (b[0], b[1], b[2], b[3], b[4], b[5]))
    outfile.close()


if __name__ == '__main__':
    civic1()
    civic2()
    civic3()
    civic4()
    civic5()
    civic6()
    civic7()
    civic8()
    civic9()

