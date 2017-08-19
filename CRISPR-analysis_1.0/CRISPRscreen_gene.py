# this script is used to calculate the statisitics of CRISPR screen at gene-level
# input sgRNA level statistics output flat file that reports statistics at the gene level

# input file:
# pickle dumped dictionary
# {gene1:{condition1:{sgRNA1:xx, sgRNA2:xx, ..},condition2:{}, ..}, gene2:{}, ..}
# 'Normalzied_RelativeAbundanceChange' is deposited in this dictionary

# pickle dumped list
# gene list in the aforementioned dictionary

# pickle dumped list
# sgRNA list in the aforementioned dictionary

# flat file of gene sgRNA position information
# sgRNA gene ORF_position(%)
# tab deliminated

# we can choose either conventional or optimized approaches for hit gene calling:
# all: all sgRNAs will be used to calculate the statistical significance
# position: search for the statistical significance peak with at least 5 sgRNAs based on subsampling based on position priority.

import os
import sys
import pickle
import numpy as np
import scipy
import math
from random import shuffle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import interp1d
from scipy.stats import mannwhitneyu
import warnings
warnings.simplefilter("error")
# alternative in mannwhitney
Testapproach='two-sided'

geneDic=sys.argv[1]
geneZDic=sys.argv[2]
geneLst=sys.argv[3]
sgRNALst=sys.argv[4]
stressed_conditionLst=sys.argv[5]
gene_sgRNA_position=sys.argv[6]
# method to set the control sgRNA set: 'all or NC'
control_setting=sys.argv[7]
if control_setting not in ['all','NC']:
    print 'incorrect control setting!'
    sys.exit()
# method for hit gene calling
hit_gene_calling=sys.argv[8]
if hit_gene_calling not in ['all','position']:
    print 'incorrect hit gene calling method setting!! please type in all or position'
    sys.exit()
FDR_threshold=float(sys.argv[9])
prefix=sys.argv[10]
sgRNAperQuasiGene=sys.argv[11]
if len(sgRNAperQuasiGene.split(','))!=2:
    print('please input two integers')
    sys.exit()
# minimal number os sgRNAs for one quasi gene
sgRNAperQuasiDown=int(sgRNAperQuasiGene.split(',')[0])
# maxium number os sgRNAs for one quasi gene
sgRNAperQuasiUp=int(sgRNAperQuasiGene.split(',')[1])
# quasi gene number for one condition
Quasi_number=int(sys.argv[12])
# minimal sgRNA number for p Value calculation in subsampling
min_sgRNAnumber=int(sys.argv[13])

Processed_geneDic=pickle.load(open(geneDic,'rb'))
Processed_geneZDic=pickle.load(open(geneZDic,'rb'))
processed_geneLst=pickle.load(open(geneLst,'rb'))
processed_sgRNALst=pickle.load(open(sgRNALst,'rb'))
stressed_conditionLst=pickle.load(open(stressed_conditionLst,'rb'))

# /////////////////////////////////////////////////////////////////
# process the gene sgRNA position flat file to construct a dictionary
# {gene:{sgRNA:position, sgRNA:position, ..}, gene:{}, ...}
gene_sgRNApositionDic={}

# used to store the sgRNA list by position order
# {gene:[sgRNA,sgR A, ..],gene:[], ..}
gene_sgRNApositionOrderDic={}
for gene in processed_geneLst:
    gene_sgRNApositionDic[gene]={}
    gene_sgRNApositionOrderDic[gene]={}
f=open(gene_sgRNA_position,'r')
for line in f:
    row=line.rstrip().split('\t')
    gene=row[0]
    sgRNA=row[1]
    position=float(row[2])
    if gene in processed_geneLst and sgRNA in processed_sgRNALst:
        gene_sgRNApositionDic[gene][sgRNA]=position
f.close()

# sort the sgRNA by position and store into gene_sgRNApositionOrderDic
for gene in processed_geneLst:
    gene_sgRNApositionOrderDic[gene]=sorted(gene_sgRNApositionDic[gene].keys(), key=lambda x:gene_sgRNApositionDic[gene][x], reverse=False)

# //////////////////////////////////////////////////////////////////
# based on the negative control sgRNAs, calculate and fit the FDR curve (function) for each condition
# to extract sub from dic based on given list
def DicExtract(sgRNALst,sgRNA_dic):
    itemLst=[]
    for sgRNA in sgRNALst:
        item=sgRNA_dic[sgRNA]
        itemLst.append(item)
    return itemLst

def quasiGene_construct(NCsgRNADic,sgRNAperQuasiDown,sgRNAperQuasiUp,Quasi_number,Testapproach):
    FDRdic={}
    NCsgRNALst=[]
    NC_relative_abundance_changeLst=[]
    for sgRNA in NCsgRNADic:
        NCsgRNALst.append(sgRNA)
        NC_relative_abundance_changeLst.append(NCsgRNADic[sgRNA])
    for sgRNAnumber in range(sgRNAperQuasiDown,(sgRNAperQuasiUp+1)):
        FDRdic[sgRNAnumber]=[]
        # ////////////// calculate random subsampling based parameters
        # subsampling for numberCross times
        for gene in range(Quasi_number):
            # creat a quasi gene 
            shuffle(NCsgRNALst)
            sgRNALst_shuffle=NCsgRNALst[:sgRNAnumber]
            #print sgRNALst_shuffle
            ChangeLst_shuffle=DicExtract(sgRNALst_shuffle,NCsgRNADic)
            # calculate the parameter
            Ustatistic,pValue_shuffle=mannwhitneyu(ChangeLst_shuffle,NC_relative_abundance_changeLst,True,Testapproach)
            Log10_pValue=math.log(pValue_shuffle,10)
            FDRdic[sgRNAnumber].append(Log10_pValue)
        FDRdic[sgRNAnumber].sort()
    return FDRdic

# use the quasiGene_construct function to construct quasiGene pValue list for each condition
# to store the MWU test Pvalue list of quasi genes
NCsgRNA_pValueDic={}
# {condition1:{sgRNAnumber:[Pvalue,pValue, ..], sgRNAnumber:[], ..}, condition2:{}, ..}
NCsgRNA_RelativeAbundanceChangeDic={}
# {condition1:[Change1,Change2, ..], condition2:[], ..}
for condition in stressed_conditionLst:
    NCsgRNA_pValueDic[condition]={}
    NCsgRNADic={}
    NCsgRNA_RelativeAbundanceChangeDic[condition]=[]
    for gene in Processed_geneDic:
        for sgRNA in Processed_geneDic[gene][condition]:
            if control_setting=='NC': # use negative control sgRNA set for normalization
                if gene=='0':
                    relative_abundance_change=Processed_geneDic[gene][condition][sgRNA]
                    NCsgRNADic[sgRNA]=relative_abundance_change
                    NCsgRNA_RelativeAbundanceChangeDic[condition].append(relative_abundance_change)
            elif control_setting=='all':
                relative_abundance_change=Processed_geneDic[gene][condition][sgRNA]
                NCsgRNADic[sgRNA]=relative_abundance_change
                NCsgRNA_RelativeAbundanceChangeDic[condition].append(relative_abundance_change)
    FDRdic=quasiGene_construct(NCsgRNADic,sgRNAperQuasiDown,sgRNAperQuasiUp,Quasi_number,Testapproach)
    # copy the result to the dictionary
    for sgRNAnumber in FDRdic:
        tempLst=FDRdic[sgRNAnumber]
        tempLst.append(-30)
        tempLst.append(0)
        NCsgRNA_pValueDic[condition][sgRNAnumber]=tempLst

# Process the NCsgRNA_pValueDic to calculate the FDR fitting and do plotting (CDF)
NCsgRNA_FDRDic={}
# {condition1:{sgRNAnumber:{'Interpolation':f,'Sorted_pValueLst':[],'CumulationLst':[]},Condition2:{}, ..}}
# Interpolation is a function used to produce FDR given a Log10pValue by interpolation 1d
cm_subsection = np.linspace(0, 0.4, 8)
colors = [ cm.jet(x) for x in cm_subsection ] 
for condition in stressed_conditionLst:
    NCsgRNA_FDRDic[condition]={}
    sgRNAnumber_strLst=[]
    for sgRNAnumber in NCsgRNA_pValueDic[condition]:
        NCsgRNA_FDRDic[condition][sgRNAnumber]={}
        Xaxis=np.sort(NCsgRNA_pValueDic[condition][sgRNAnumber])
        NCsgRNA_FDRDic[condition][sgRNAnumber]['Sorted pValueLst']=Xaxis
        N=len(NCsgRNA_pValueDic[condition][sgRNAnumber])
        if N!=(Quasi_number+2):
            print 'not enough quasi genes!'
            sys.exit()
        else:
            Yaxis = np.array(range(N))/(float(N)-1)
            NCsgRNA_FDRDic[condition][sgRNAnumber]['CumulationLst']=Yaxis
            interpolationFunction = interp1d(Xaxis, Yaxis)
            NCsgRNA_FDRDic[condition][sgRNAnumber]['Interpolation']=interpolationFunction
        if sgRNAnumber in [3,5,10,15,20,25,30,40]:
            sgRNAnumber_strLst.append('sgRNAnumber= '+str(sgRNAnumber))
            color=colors[[3,5,10,15,20,25,30,40].index(sgRNAnumber)]
            plt.plot(Xaxis, Yaxis,color=color)
    plt.axhline(y=FDR_threshold,xmin=-5,xmax=1,linewidth=2, color='#000000',linestyle='--')
    sgRNAnumber_strLst.append('FDR threshold= %s'%(FDR_threshold))
    plt.legend(sgRNAnumber_strLst, loc='best')
    plt.xlabel('Log10 MWU test P value',fontsize=20)
    plt.ylabel('False discovery rate',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((-3,-0.5))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.95)
    plt.savefig('%s_results/%s_%s_quasigeneFDR.png'%(prefix,prefix,condition),dpi=400)
    plt.clf()
os.system('mkdir %s_results/%s_quasigeneFDR/'%(prefix,prefix))
os.system('mv %s_results/*_quasigeneFDR.png %s_results/%s_quasigeneFDR/'%(prefix,prefix,prefix))

# ///////////////////////////////////////////////////////////////////////////
# start to process the gene level statistics calculation
Result_Dic={}
# {gene:{'MedianRAC':xx,'MeidanZ':xx,'Log10Pvalue':xx,'sgRNAnumber':xx,'FDRvalue':xx}, gene:{}, ..}


for condition in stressed_conditionLst:
    Result_Dic[condition]={}
    for gene in processed_geneLst:
        if gene!='0':
            Result_Dic[condition][gene]={}
            NC_relativechangeLst=NCsgRNA_RelativeAbundanceChangeDic[condition]
            if hit_gene_calling=='all':
                sgRNAnumber=len(Processed_geneDic[gene][condition])
                sgRNALst=gene_sgRNApositionOrderDic[gene][:sgRNAnumber]
                relativechangeLst=DicExtract(sgRNALst,Processed_geneDic[gene][condition])
                ZscoreLst=DicExtract(sgRNALst,Processed_geneZDic[gene][condition])
                medianRAC=np.median(np.array(relativechangeLst))
                medianZ=np.median(np.array(ZscoreLst))
                Ustatistic,pValue=mannwhitneyu(relativechangeLst,NC_relativechangeLst,True,Testapproach)
                Log10P=math.log(pValue,10)
                FDRvalue=NCsgRNA_FDRDic[condition][sgRNAnumber]['Interpolation'](Log10P)
                Result_Dic[condition][gene]['MedianRAC']=medianRAC
                Result_Dic[condition][gene]['MedianZ']=medianZ
                Result_Dic[condition][gene]['Log10Pvalue']=Log10P
                Result_Dic[condition][gene]['sgRNAnumber']=sgRNAnumber
                Result_Dic[condition][gene]['FDRvalue']=FDRvalue
            # if position-priority based subsampling is used, then search the peak of significance and store the values into the dic if better results detected
            elif hit_gene_calling=='position':
                available_sgRNAnumber=min(len(Processed_geneDic[gene][condition]),sgRNAperQuasiUp)
                Result_Dic[condition][gene]['MedianRAC']=0
                Result_Dic[condition][gene]['MedianZ']=0
                Result_Dic[condition][gene]['Log10Pvalue']=1
                Result_Dic[condition][gene]['sgRNAnumber']=sgRNAperQuasiDown
                Result_Dic[condition][gene]['FDRvalue']=1
                if min_sgRNAnumber<available_sgRNAnumber+1:
                    temporary_sgRNAnumber=min_sgRNAnumber
                else:
                    print('The gene %s does not has %s sgRNAs'%(gene,min_sgRNAnumber))
                    temporary_sgRNAnumber=1
                for sgRNAnumber in range(temporary_sgRNAnumber,available_sgRNAnumber+1):
                    sgRNALst_position=gene_sgRNApositionOrderDic[gene][:sgRNAnumber]
                    relativechangeLst=DicExtract(sgRNALst_position,Processed_geneDic[gene][condition])
                    ZscoreLst=DicExtract(sgRNALst_position,Processed_geneZDic[gene][condition])
                    medianRAC=np.median(np.array(relativechangeLst))
                    medianZ=np.median(np.array(ZscoreLst))
                    Ustatistic,pValue=mannwhitneyu(relativechangeLst,NC_relativechangeLst,True,Testapproach)
                    Log10P=math.log(pValue,10)
                    FDRvalue=NCsgRNA_FDRDic[condition][sgRNAnumber]['Interpolation'](Log10P)
                    if Log10P<Result_Dic[condition][gene]['Log10Pvalue']:
                        Result_Dic[condition][gene]['MedianRAC']=medianRAC
                        Result_Dic[condition][gene]['MedianZ']=medianZ
                        Result_Dic[condition][gene]['Log10Pvalue']=Log10P
                        Result_Dic[condition][gene]['sgRNAnumber']=sgRNAnumber
                        Result_Dic[condition][gene]['FDRvalue']=FDRvalue

pickle.dump(Result_Dic,open('%s_results/%s_Result_geneDic.pickle'%(prefix,prefix),'wb'))  #


# ///////////////////////////////////////////////////////////////////
# output all results into a csv file
os.system('mkdir %s_results/%s_gene_statistics/'%(prefix,prefix))
for condition in stressed_conditionLst:
    os.system('cat /dev/null > %s_results/%s_gene_statistics/%s_%s_gene_statistics.txt'%(prefix,prefix,prefix,condition))
    g=open('%s_results/%s_gene_statistics/%s_%s_gene_statistics.txt'%(prefix,prefix,prefix,condition),'r+')
    g.write('gene\tMedianRAC\tMedianZ\t-Log10Pvalue\tsgRNAnumber\tFDRvalue\n')
    for gene in processed_geneLst:
        if gene!='0':
            medianRAC=Result_Dic[condition][gene]['MedianRAC']
            medianZ=Result_Dic[condition][gene]['MedianZ']
            Log10P=Result_Dic[condition][gene]['Log10Pvalue']
            sgRNAnumber=Result_Dic[condition][gene]['sgRNAnumber']
            FDRvalue=Result_Dic[condition][gene]['FDRvalue']
            g.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(gene,str(medianRAC),str(medianZ),str(-1*Log10P),str(sgRNAnumber),str(FDRvalue)))
    g.close()
