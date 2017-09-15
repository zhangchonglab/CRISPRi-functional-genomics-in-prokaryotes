#Operon file
#operon1 name1 gene
#operon2 name2 gene .. .. ..

import os
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt

Finish_result_geneDic=sys.argv[1]
Operon_gene_List=sys.argv[2]
prefix=sys.argv[3]
Result_geneDic=pickle.load(open(Finish_result_geneDic,'rb'))

# used to store the Operon list by gene
#{Operon1:[gene1,gene2],Operon2:[gene3,gene4]}
OperonDic={}  #{operon1:[gene1 ,gene2],operon2:[gene0,gene1],...}

#extract the information of Operon 
with open(Operon_gene_List,'r') as f:
    for (i,line) in enumerate(f):
        row=line.rstrip().split('\t')
        if i!=0:
            OperonDic[row[0]]=row[1].split(',')

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
for operon in OperonDic:
    for condition in Result_geneDic:
        Result_OperonDic[condition][operon]={}
        for gene in Result_geneDic[condition]:
            if gene in OperonDic[operon]:
                Result_OperonDic[condition][operon][gene]=[]
                medianRAC=Result_geneDic[condition][gene]['MedianRAC']
                medianZ=Result_geneDic[condition][gene]['MedianZ']
                Log10P=Result_geneDic[condition][gene]['Log10Pvalue']
                sgRNAnumber=Result_geneDic[condition][gene]['sgRNAnumber']
                FDRvalue=Result_geneDic[condition][gene]['FDRvalue']
                Result_OperonDic[condition][operon][gene].append(medianRAC)
                Result_OperonDic[condition][operon][gene].append(medianZ)
                Result_OperonDic[condition][operon][gene].append(Log10P)
                Result_OperonDic[condition][operon][gene].append(sgRNAnumber)
                Result_OperonDic[condition][operon][gene].append(FDRvalue)
                Result_geneLst[condition].append(gene)

#store the gene that is not in operon
for condition in Result_OperonDic:
    for gene in Result_geneDic[condition]:
        if gene not in Result_geneLst[condition]:
            RemoveGene[condition].append(gene)

#remove the empty operon
for condition in Result_OperonDic:
    for operon in OperonDic:
        if len(Result_OperonDic[condition][operon])==0:
            del Result_OperonDic[condition][operon]
            RemoveOperon[condition].append(operon)

# /////////////////////////////////
#output all results into a file
# ////////////////////////////////
os.system('mkdir %s_results/%s_operon_statistics/'%(prefix,prefix))
#output the operon structure
for condition in Result_OperonDic:
    os.system('cat /dev/null > %s_results/%s_operon_statistics/%s_%s_operon_statistics.txt'%(prefix,prefix,prefix,condition))
    with open('%s_results/%s_operon_statistics/%s_%s_operon_statistics.txt'%(prefix,prefix,prefix,condition),'r+') as f: 
        f.write('%s\nGene\tsgRNAnumber\tmedianRAC\tmedianZ\t-Log10P\tFDRvalue\n'%(condition))
        for operon in Result_OperonDic[condition]:
            geneName=''
            sgRNAnumber=''
            medianRAC=''
            medianZ=''
            negativeLog10P=''
            FDRvalue=''
            for gene in Result_OperonDic[condition][operon]:
                geneName+=gene+','
                medianRAC+=str(Result_OperonDic[condition][operon][gene][0])+','
                medianZ+=str(Result_OperonDic[condition][operon][gene][1])+','
                negativeLog10P+=str(-1*Result_OperonDic[condition][operon][gene][2])+','
                sgRNAnumber+=str(Result_OperonDic[condition][operon][gene][3])+','
                FDRvalue+=str(Result_OperonDic[condition][operon][gene][4])+','
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(geneName.strip(','),sgRNAnumber.strip(','),medianRAC.strip(','),medianZ.strip(','),negativeLog10P.strip(','),FDRvalue.strip(',')))

#output the RemoveGene
for condition in RemoveGene:
    if len(RemoveGene[condition])!=0:
        os.system('cat /dev/null > %s_results/%s_operon_statistics/%s_%s_RemoveGene_statistics.txt' %(prefix,prefix,prefix,condition))
        with open('%s_results/%s_operon_statistics/%s_%s_RemoveGene_statistics.txt' %(prefix,prefix,prefix,condition),'r+') as f:
            f.write('%s\nGene'%(condition))
            for line in RemoveGene[condition]:
                f.write('%s\n'%(line))

#output the RemoveOperon
for condition in RemoveOperon:
    if len(RemoveOperon[condition])!=0:
        os.system('cat /dev/null > %s_results/%s_operon_statistics/%s_%s_RemoveOperon_statistics.txt' %(prefix,prefix,prefix,condition))
        with open('%s_results/%s_operon_statistics/%s_%s_RemoveOperon_statistics.txt' %(prefix,prefix,prefix,condition),'r+') as f:
            f.write('%s\nOperon'%(condition))
            for line in RemoveOperon[condition]:
                f.write('%s\n'%(line))
