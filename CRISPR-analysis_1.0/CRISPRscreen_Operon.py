

#Operon file
#operon1 name1 gene
#operon2 name2 gene .. .. ..

import os
import sys
import pickle
import numpy as np
import scipy
import math
from random import shuffle
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import interp1d
from scipy.stats import mannwhitneyu
import warnings
warnings.simplefilter("error")

Finish_result_geneDic=sys.argv[1]
Operon_gene_List=sys.argv[2]
prefix=sys.argv[3]
Result_geneDic=pickle.load(open(Finish_result_geneDic,'rb'))


# used to store the Operon list by gene
#{Operon1:[gene1,gene2],Operon2:[gene3,gene4]}
OperonLst={}  #{operon1:operon1 name,operon2:operon2 name,...}
OperonDic={}  #{operon1:[gene1 ,gene2],operon2:[gene0,gene1],...}

#extract the information of Operon 
with open(Operon_gene_List,'r') as f:
    for (i,line) in enumerate(f):
        row=line.rstrip().split('\t')
        if i!=0:
            OperonLst[row[0]]=row[1]
            OperonDic[row[0]]=[]
            lenth=len(row[2].split(','))
            for i in range(lenth):
                OperonDic[row[0]].append(row[2].split(',')[i])


#Integration the result gene based on the Operon
Result_geneLst={}
RemoveGene={}
RemoveOperon={}
Result_OperonDic={}
for condition in Result_geneDic:
    RemoveOperon[condition]=[]
    Result_geneLst[condition]=[]
    RemoveGene[condition]=[]
    Result_OperonDic[condition]={}
for record in OperonLst:
    for condition in Result_geneDic:
        Result_OperonDic[condition][record]={}
        for gene in Result_geneDic[condition]:
            if gene in OperonDic[record]:
                Result_OperonDic[condition][record][gene]=[]
                medianZ=Result_geneDic[condition][gene]['Median']
                Log10P=Result_geneDic[condition][gene]['Log10Pvalue']
                sgRNAnumber=Result_geneDic[condition][gene]['sgRNAnumber']
                FDRvalue=Result_geneDic[condition][gene]['FDRvalue']
                Result_OperonDic[condition][record][gene].append(medianZ)
                Result_OperonDic[condition][record][gene].append(Log10P)
                Result_OperonDic[condition][record][gene].append(sgRNAnumber)
                Result_OperonDic[condition][record][gene].append(FDRvalue)
                Result_geneLst[condition].append(gene)
#store the gene that is not in operon
for condition in Result_OperonDic:
    for line in Result_geneDic[condition]:
        if line not in Result_geneLst[condition]:
            RemoveGene[condition].append(line)

#remove the empty operon which there is no gene in 
for condition in Result_OperonDic:
    for record in OperonLst:
        if len(Result_OperonDic[condition][record])==0:
            del Result_OperonDic[condition][record]
            RemoveOperon[condition].append(record)

# /////////////////////////////////
#output all results into a file
# ////////////////////////////////
os.system('mkdir %s_results/%s_operon_statistics/'%(prefix,prefix))
#output the operon struction
for condition in Result_OperonDic:
    os.system('cat /dev/null > %s_results/%s_operon_statistics/%s_%s_operon_statistics.txt'%(prefix,prefix,prefix,condition))
    with open('%s_results/%s_operon_statistics/%s_%s_operon_statistics.txt'%(prefix,prefix,prefix,condition),'r+') as f: 
        f.write('%s\nOperon\tGene\tsgRNAnumber\tmedianZ\tLog10P\tFDRvalue\n'%(condition))
        for operon in Result_OperonDic[condition]:
            geneName=operon+'\t'
            sgRNAnumber=''
            medianZ=''
            Log10P=''
            FDRvalue=''
            for rec in Result_OperonDic[condition][operon]:
                geneName+=rec+','
                medianZ+=str(Result_OperonDic[condition][operon][rec][0])+','
                Log10P+=str(Result_OperonDic[condition][operon][rec][1])+','
                sgRNAnumber+=str(Result_OperonDic[condition][operon][rec][2])+','
                FDRvalue+=str(Result_OperonDic[condition][operon][rec][3])+','
            f.write('%s\t%s\t%s\t%s\t%s\n'%(geneName.strip(','),sgRNAnumber.strip(','),medianZ.strip(','),Log10P.strip(','),FDRvalue.strip(',')))

#output the RemoveGene
for condition in RemoveGene:
    if len(RemoveGene[condition])!=0:
        os.system('cat/dev/null > %s_results/%s_operon_statistics/%s_%s_RemoveGene_statistics.txt' %(prefix,prefix,prefix,condition))
        with open('%s_results/%s_operon_statistics/%s_%s_RemoveGene_statistics.txt' %(prefix,prefix,prefix,condition),'r+') as f:
            f.write('%s\nGene'%(condition))
            for line in RemoveGene[condition]:
                f.write('%s\n'%(line))

#output the RemoveOperon
for condition in RemoveOperon:
    if len(RemoveOperon[condition])!=0:
        os.system('cat/dev/null > %s_results/%s_operon_statistics/%s_%s_RemoveOperon_statistics.txt' %(prefix,prefix,prefix,condition))
        with open('%s_results/%s_operon_statistics/%s_%s_RemoveOperon_statistics.txt' %(prefix,prefix,prefix,condition),'r+') as f:
            f.write('%s\nOperon'%(condition))
            for line in RemoveOperon[condition]:
                f.write('%s\n'%(line))



        


