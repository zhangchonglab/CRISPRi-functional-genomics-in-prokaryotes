import os
import sys
import pickle
import numpy as np
import scipy
import math
import pandas as pd
from random import shuffle
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
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
# minimal sgRNA number for p Value calculation in subsampling
min_sgRNAnumber=int(sys.argv[12])

Processed_geneDic=pickle.load(open(geneDic,'rb'))
Processed_geneZDic=pickle.load(open(geneZDic,'rb'))
processed_geneLst=pickle.load(open(geneLst,'rb'))
processed_sgRNALst=pickle.load(open(sgRNALst,'rb'))
stressed_conditionLst=pickle.load(open(stressed_conditionLst,'rb'))
# quasi gene number for one condition
Quasi_number=min(10*len(processed_geneLst),10000)

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
# based on the negative control sgRNAs, calculate and fit the FPR curve (function) for each condition
# to extract sub from dic based on given list
def DicExtract(sgRNALst,sgRNA_dic):
    itemLst=[]
    for sgRNA in sgRNALst:
        item=sgRNA_dic[sgRNA]
        itemLst.append(item)
    return itemLst

def quasiGene_construct(NCsgRNADic,sgRNAperQuasiDown,sgRNAperQuasiUp,Quasi_number,Testapproach):
    FPRdic={}
    NCsgRNALst=[]
    NC_relative_abundance_changeLst=[]
    for sgRNA in NCsgRNADic:
        NCsgRNALst.append(sgRNA)
        NC_relative_abundance_changeLst.append(NCsgRNADic[sgRNA])
    for sgRNAnumber in range(sgRNAperQuasiDown,(sgRNAperQuasiUp+1)):
        FPRdic[sgRNAnumber]=[]
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
            negativeLog10_pValue=-math.log(pValue_shuffle,10)
            medianRAC=np.median(np.array(ChangeLst_shuffle))
            score=negativeLog10_pValue*abs(medianRAC)
            FPRdic[sgRNAnumber].append(score)
        FPRdic[sgRNAnumber].sort()
    return FPRdic

import scipy as sp
from scipy import interpolate


def Qvalue_estimate(pv, m=None, pi0=None):
    """
    Estimates q-values from p-values
    Args
    =====
    m: number of tests. If not specified m = pv.size
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1
    """
    if np.min(pv) <= 0 or np.max(pv) > 1:
        print ("p-values should be within (0,1]")
        sys.exit(-1)
    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()
    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0
    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = scipy.arange(0, 0.95, 0.01)
        counts = scipy.array([(pv > i).sum() for i in scipy.arange(0, 0.95, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))
        pi0 = scipy.array(pi0)
        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if pi0 > 1:
            pi0 = 1.0
    if pi0<0.9:
        pi0=0.9
    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0
    p_ordered = scipy.argsort(pv)
    pv = pv[p_ordered]
    qv = pi0 * m/len(pv) * pv
    qv[-1] = min(qv[-1], 1.0)
    for i in xrange(len(pv)-2, -1, -1):
        qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])
    # reorder qvalues
    qv_temp = qv.copy()
    qv = scipy.zeros_like(qv)
    qv[p_ordered] = qv_temp
    # reshape qvalues
    qv = qv.reshape(original_shape)
    return pi0,qv

def FDR_estimate(scoreLst,NCscoreLst,score):
    """
    Estimates FDR from two score lists genes and quasi genes (constructed from NC sgRNAs), respectively
    """
    sigLst=[x for x in scoreLst if x >= score]
    NC_sigLst=[x for x in NCscoreLst if x >= score]
    FDR=float(len(NC_sigLst))/( float(len(NC_sigLst)) + float(len(sigLst)))
    return FDR

# use the quasiGene_construct function to construct quasiGene pValue list for each condition
# to store the MWU test Pvalue list of quasi genes
NCsgRNA_scoreDic={}
# {condition1:{sgRNAnumber:[Pvalue,pValue, ..], sgRNAnumber:[], ..}, condition2:{}, ..}
NCsgRNA_RelativeAbundanceChangeDic={}
# {condition1:[Change1,Change2, ..], condition2:[], ..}
for condition in stressed_conditionLst:
    NCsgRNA_scoreDic[condition]={}
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
    FPRdic=quasiGene_construct(NCsgRNADic,sgRNAperQuasiDown,sgRNAperQuasiUp,Quasi_number,Testapproach)
    # copy the result to the dictionary
    for sgRNAnumber in FPRdic:
        tempLst=FPRdic[sgRNAnumber]
        tempLst.append(0.0)
        tempLst.append(200.0)
        NCsgRNA_scoreDic[condition][sgRNAnumber]=tempLst

# Process the NCsgRNA_scoreDic to calculate the FPR fitting and do plotting (CDF)
# store score-FPR function
NCsgRNA_FPRDic={}
# store FPR-score function
NCsgRNA_FPRrevDic={}
# {condition1:{sgRNAnumber:{'Interpolation':f,'Sorted scoreLst':[],}Condition2:{}, ..}}
# Interpolation is a function used to produce FPR given a Log10pValue by interpolation 1d
cm_subsection = np.linspace(0, 0.4, 10)
colors = [ cm.jet(x) for x in cm_subsection ] 
for condition in stressed_conditionLst:
    NCsgRNA_FPRDic[condition]={}
    NCsgRNA_FPRrevDic[condition]={}
    sgRNAnumber_strLst=[]
    for sgRNAnumber in NCsgRNA_scoreDic[condition]:
        NCsgRNA_FPRDic[condition][sgRNAnumber]={}
        NCsgRNA_FPRrevDic[condition][sgRNAnumber]={}
        Xaxis=np.sort(NCsgRNA_scoreDic[condition][sgRNAnumber])
        NCsgRNA_FPRDic[condition][sgRNAnumber]['Sorted scoreLst']=Xaxis
        N=len(NCsgRNA_scoreDic[condition][sgRNAnumber])
        if N!=(Quasi_number+2):
            print ('not enough quasi genes!')
            sys.exit()
        else:
            Yaxis = 1-np.array(range(N))/(float(N)-1)
            interpolationFunction = interp1d(Xaxis, Yaxis)
            Rev_interpolationFunction = interp1d(Yaxis, Xaxis)
            NCsgRNA_FPRDic[condition][sgRNAnumber]['Interpolation']=interpolationFunction
            NCsgRNA_FPRrevDic[condition][sgRNAnumber]['Interpolation']=Rev_interpolationFunction
        if sgRNAnumber in [1,2,3,5,10,15,20]:
            sgRNAnumber_strLst.append('sgRNAnumber= '+str(sgRNAnumber))
            color=colors[[1,2,3,5,10,15,20].index(sgRNAnumber)]
            plt.plot(Xaxis, Yaxis,color=color)
    plt.legend(sgRNAnumber_strLst, loc='best')
    plt.xlabel('$-Log_{10}$ MWU test $P * abs(medianRAC)$ value',fontsize=14)
    plt.ylabel('False positive rate (FPR)',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((0,5))
    plt.ylim((0,0.3))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    plt.savefig('%s_results/%s_%s_quasigeneFPR.png'%(prefix,prefix,condition),dpi=400)
    plt.clf()
os.system('mkdir -p %s_results/%s_quasigeneFPR/'%(prefix,prefix))
os.system('mv %s_results/*_quasigeneFPR.png %s_results/%s_quasigeneFPR/'%(prefix,prefix,prefix))

# store the score cutoff for several commonly used FPR
for sgRNAnumber in range(sgRNAperQuasiDown,(sgRNAperQuasiUp+1)):
    myFPRbook={}
    for condition in stressed_conditionLst:
        myFPRbook[condition]={}
        for FPR in [0.001,0.005,0.01,0.02,0.05,0.1]:
            score=NCsgRNA_FPRrevDic[condition][sgRNAnumber]['Interpolation'](FPR)
            myFPRbook[condition][FPR]=score
    myFPR_DF=pd.DataFrame(myFPRbook)
    myFPR_DF['FPR']=myFPR_DF.index
    myFPR_DF.to_csv('%s_results/%s_quasigeneFPR/FPR_sgRNAnumber=%d.txt'%(prefix,prefix,sgRNAnumber),index=None,sep='\t')

# ///////////////////////////////////////////////////////////////////////////
# start to process the gene level statistics calculation
# {gene:{'MedianRAC':xx,'MeidanZ':xx,'Log10PvalueMWU':xx,'Log10PvalueT':xx,'sgRNAnumber':xx,'score':xx, 'FDRvalue':xx, 'Qvalue':xx, }, gene:{}, ..}
Result_Dic={}
# record all T test P value for genes
gene_Ttest_pValue_Dic={}
# record all FPR value for genes
gene_FPR_Dic={}
# record scores of all genes: {condition:[], condition:[]}
gene_score_Dic={}
# record scores of all quasi genes: {condition:[], condition:[]}
NC_score_Dic={}

# ============================================
# Calculate T test, MWU test for each gene; construct quasi gene according to gene profile (sgRNA number); record relevant values
for condition in stressed_conditionLst:
    Result_Dic[condition]={}
    gene_Ttest_pValue_Dic[condition]=[]
    gene_FPR_Dic[condition]=[]
    gene_score_Dic[condition]=[]
    NC_score_Dic[condition]=[]
    for gene in processed_geneLst:
        if gene!='0':
            Result_Dic[condition][gene]={}
            NC_relativechangeLst=NCsgRNA_RelativeAbundanceChangeDic[condition]
            if hit_gene_calling=='all':
                sgRNAnumber=len(Processed_geneDic[gene][condition])
                if sgRNAnumber>=1:
                    # get all available sgRNAs for this gene passing the control condition read number threshold
                    sgRNALst=[x for x in Processed_geneDic[gene][condition]]
                    relativechangeLst=DicExtract(sgRNALst,Processed_geneDic[gene][condition])
                    ZscoreLst=DicExtract(sgRNALst,Processed_geneZDic[gene][condition])
                    medianRAC=np.median(np.array(relativechangeLst))
                    medianZ=np.median(np.array(ZscoreLst))
                    Ustatistic,pValue=mannwhitneyu(relativechangeLst,NC_relativechangeLst,True,Testapproach)
                    Log10P_MWUtest=math.log(pValue,10)
                    score=abs(medianRAC)*(-1*Log10P_MWUtest)
                    FPRvalue=NCsgRNA_FPRDic[condition][sgRNAnumber]['Interpolation'](score)
                    if sgRNAnumber>=2:
                        Tstat,Ttest_pValue=ttest_ind(np.array(relativechangeLst), np.array(NC_relativechangeLst))
                    else: # in the case only one sgRNA is available, replicate it to perform the T test
                        relativechange=relativechangeLst[0]
                        relativechangeLst.append(relativechange)
                        Tstat,Ttest_pValue=ttest_ind(np.array(relativechangeLst), np.array(NC_relativechangeLst))
                else:
                    medianRAC=0; medianZ=0; Log10P_MWUtest=0; Ttest_pValue=1; sgRNAnumber=0; FPRvalue=1
                Result_Dic[condition][gene]['MedianRAC']=medianRAC
                Result_Dic[condition][gene]['MedianZ']=medianZ
                Result_Dic[condition][gene]['Log10PvalueMWU']=Log10P_MWUtest
                Result_Dic[condition][gene]['PvalueTtest']=Ttest_pValue
                Result_Dic[condition][gene]['sgRNAnumber']=sgRNAnumber
                Result_Dic[condition][gene]['score']=score
                Result_Dic[condition][gene]['FPRvalue']=FPRvalue
            # if position-priority based subsampling is used, then search the peak of significance and store the values into the dic if better results detected
            elif hit_gene_calling=='position':
                available_sgRNAnumber=min(len(Processed_geneDic[gene][condition]),sgRNAperQuasiUp)
                Result_Dic[condition][gene]['MedianRAC']=0
                Result_Dic[condition][gene]['MedianZ']=0
                Result_Dic[condition][gene]['Log10PvalueMWU']=0
                Result_Dic[condition][gene]['PvalueTtest']=1
                Result_Dic[condition][gene]['sgRNAnumber']=0
                Result_Dic[condition][gene]['score']=0
                Result_Dic[condition][gene]['FPRvalue']=1
                if available_sgRNAnumber>=1:
                    # get all available sgRNAs for this gene passing the control condition read number threshold
                    mysgRNALst=[x for x in Processed_geneDic[gene][condition]]
                    # get the sgRNAs from the position ordered list and filter those <=readthreshold in control condition
                    thisgene_sgRNApositionOrderLst=[sgRNA for sgRNA in gene_sgRNApositionOrderDic[gene] if sgRNA in mysgRNALst]
                    if min_sgRNAnumber<available_sgRNAnumber+1:
                        temporary_sgRNAnumber=min_sgRNAnumber
                    else:
                        temporary_sgRNAnumber=1
                    for sgRNAnumber in range(temporary_sgRNAnumber,available_sgRNAnumber+1):
                        sgRNALst_position=thisgene_sgRNApositionOrderLst[:sgRNAnumber]
                        relativechangeLst=DicExtract(sgRNALst_position,Processed_geneDic[gene][condition])
                        ZscoreLst=DicExtract(sgRNALst_position,Processed_geneZDic[gene][condition])
                        medianRAC=np.median(np.array(relativechangeLst))
                        medianZ=np.median(np.array(ZscoreLst))
                        Ustatistic,pValue=mannwhitneyu(relativechangeLst,NC_relativechangeLst,True,Testapproach)
                        Log10P_MWUtest=math.log(pValue,10)
                        score=abs(medianRAC)*(-1*Log10P_MWUtest)
                        FPRvalue=NCsgRNA_FPRDic[condition][sgRNAnumber]['Interpolation'](score)
                        if sgRNAnumber>=2:
                            Tstat,Ttest_pValue=ttest_ind(np.array(relativechangeLst), np.array(NC_relativechangeLst))
                        else: # in the case only one sgRNA is available, replicate it to perform the T test
                            relativechange=relativechangeLst[0]
                            relativechangeLst.append(relativechange)
                            Tstat,Ttest_pValue=ttest_ind(np.array(relativechangeLst), np.array(NC_relativechangeLst))
                        if FPRvalue<Result_Dic[condition][gene]['FPRvalue']:
                            Result_Dic[condition][gene]['score']=score
                            Result_Dic[condition][gene]['MedianRAC']=medianRAC
                            Result_Dic[condition][gene]['MedianZ']=medianZ
                            Result_Dic[condition][gene]['Log10PvalueMWU']=Log10P_MWUtest
                            Result_Dic[condition][gene]['PvalueTtest']=Ttest_pValue
                            Result_Dic[condition][gene]['sgRNAnumber']=sgRNAnumber
                            Result_Dic[condition][gene]['FPRvalue']=FPRvalue
            if Result_Dic[condition][gene]['sgRNAnumber']>=1:
                sgRNAnumber=Result_Dic[condition][gene]['sgRNAnumber']; Ttest_pValue=Result_Dic[condition][gene]['PvalueTtest']; score=Result_Dic[condition][gene]['score']; FPRvalue=Result_Dic[condition][gene]['FPRvalue']
                myNC_scoreLst=NCsgRNA_scoreDic[condition][sgRNAnumber]
                shuffle(myNC_scoreLst)
                if myNC_scoreLst[0]!=200 and myNC_scoreLst[0]!=0:
                    NC_score_Dic[condition].append(myNC_scoreLst[0])
                gene_Ttest_pValue_Dic[condition].append(Ttest_pValue)
                gene_score_Dic[condition].append(score)
                gene_FPR_Dic[condition].append(FPRvalue)
# /////////////////////////////////////////////////////////
# calculate the FDR value according to quasi genes and Q value according to T test P value

# ==================================
# overall FDR evaluation by constructing quasi genes
gene_FDRDic={}
gene_FDRrevDic={}
for condition in stressed_conditionLst:
    gene_FDRDic[condition]=0
    gene_FDRrevDic[condition]=0
    score_sortLst=gene_score_Dic[condition]
    score_sortLst.sort()
    Xaxis=[0.0]; Yaxis=[1.0]
    for score in score_sortLst:
        FDR=FDR_estimate(gene_score_Dic[condition],NC_score_Dic[condition],score)
        Xaxis.append(score); Yaxis.append(FDR)
    Xaxis.append(200.0); Yaxis.append(0.0)
    interpolationFunction = interp1d(np.array(Xaxis), np.array(Yaxis))
    Rev_interpolationFunction = interp1d(np.array(Yaxis), np.array(Xaxis))
    gene_FDRDic[condition]=interpolationFunction
    gene_FDRrevDic[condition]=Rev_interpolationFunction
    plt.plot(Xaxis, Yaxis,color='#82E0AA')
    plt.axhline(y=FDR_threshold, linewidth=2, color='#000000',linestyle='--')
    plt.text(2,0.18,'FDR threshold = %s'%(FDR_threshold),fontsize=12)
    plt.xlabel('$-Log_{10}$ MWU test $P * abs(medianRAC)$ value',fontsize=14)
    plt.ylabel('False discovery rate',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((0,4))
    plt.ylim((0,0.2))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    plt.savefig('%s_results/%s_%s_quasigeneFDR.png'%(prefix,prefix,condition),dpi=400)
    plt.clf()
os.system('mkdir -p %s_results/%s_quasigeneFDR/'%(prefix,prefix))
os.system('mv %s_results/*_quasigeneFDR.png %s_results/%s_quasigeneFDR/'%(prefix,prefix,prefix))

# store the score cutoff for several commonly used FDR
myFDRbook={}
for condition in stressed_conditionLst:
    myFDRbook[condition]={}
    for FDR in [0.001,0.005,0.01,0.02,0.05,0.1]:
        score=gene_FDRrevDic[condition](FDR)
        myFDRbook[condition][FDR]=score
    myFDR_DF=pd.DataFrame(myFDRbook)
    myFDR_DF['FDR']=myFDR_DF.index
    myFDR_DF.to_csv('%s_results/%s_quasigeneFDR/FDR_score.txt'%(prefix,prefix),index=None,sep='\t')

# ==================================
# draw plot of T test P value or FPR value vs. Q value (Storey Tibshirani approach) relationship
gene_QTp_valueDic={}
gene_QTp_valuerevDic={}
gene_QFPR_valueDic={}
gene_QFPR_valuerevDic={}
bins=np.linspace(0,1,50)
for condition in stressed_conditionLst:
    gene_QTp_valueDic[condition]=0
    gene_QTp_valuerevDic[condition]=0
    gene_QFPR_valueDic[condition]=0
    gene_QFPR_valuerevDic[condition]=0
    pv_T=np.sort(np.array(gene_Ttest_pValue_Dic[condition]))
    plt.hist(pv_T,bins)
    plt.xlabel('T test $P$ value',fontsize=18)
    plt.ylabel('gene number',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    plt.savefig('%s_results/%s_%s_Ttest_pValue.png'%(prefix,prefix,condition),dpi=400)
    plt.clf()
    pv_FPR=np.sort(np.array(gene_FPR_Dic[condition]))
    plt.hist(pv_FPR,bins)
    plt.xlabel('Score based $FPR$ value',fontsize=18)
    plt.ylabel('gene number',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    plt.savefig('%s_results/%s_%s_scoreFPR.png'%(prefix,prefix,condition),dpi=400)
    plt.clf()
    # pi0 calculated here will be used in Q value calculation based on FPR later
    pi0,qv_T=Qvalue_estimate(pv_T)
    Xaxis_T=pv_T.tolist(); Yaxis_T=qv_T.tolist(); Xaxis_T.append(0.0); Yaxis_T.append(0.0); Xaxis_T.append(1.0); Yaxis_T.append(1.0)
    interpolationFunction = interp1d(np.array(Xaxis_T), np.array(Yaxis_T))
    Rev_interpolationFunction = interp1d(np.array(Yaxis_T), np.array(Xaxis_T))
    gene_QTp_valueDic[condition]=interpolationFunction
    gene_QTp_valuerevDic[condition]=Rev_interpolationFunction
    # use pi0 from T test P value derived Storey approach
    pi0new,qv_FPR=Qvalue_estimate(pv_FPR,pi0=pi0)
    Xaxis_FPR=pv_FPR.tolist(); Yaxis_FPR=qv_FPR.tolist(); Xaxis_FPR.append(0.0); Yaxis_FPR.append(0.0); Xaxis_FPR.append(1.0); Yaxis_FPR.append(1.0)
    interpolationFunction = interp1d(np.array(Xaxis_FPR), np.array(Yaxis_FPR))
    Rev_interpolationFunction = interp1d(np.array(Yaxis_FPR), np.array(Xaxis_FPR))
    gene_QFPR_valueDic[condition]=interpolationFunction
    gene_QFPR_valuerevDic[condition]=Rev_interpolationFunction
    plt.plot(pv_T,qv_T,color='#58D68D',ls='--',label='T test $P$ value')
    plt.plot(pv_FPR,qv_FPR,color='#5DADE2',ls='--',label='Score based $FPR$ value')
    plt.legend(loc='upper left',fontsize=12)
    plt.axhline(y=FDR_threshold, linewidth=2, color='#000000',linestyle='--')
    plt.text(0.001,0.15,'FDR threshold = %s'%(FDR_threshold),fontsize=12)
    plt.xlabel('$P$ value',fontsize=18)
    plt.ylabel('$Q$ value',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((0,0.1))
    plt.ylim((0,0.2))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    plt.savefig('%s_results/%s_%s_Pvalue_Qvalue.png'%(prefix,prefix,condition),dpi=400)
    plt.clf()
os.system('mkdir -p %s_results/%s_Pvalue_Qvalue/'%(prefix,prefix))
os.system('mv %s_results/*.png %s_results/%s_Pvalue_Qvalue/'%(prefix,prefix,prefix))

# store the T test P value cutoff for several commonly used Q value
myQTp_book={}
for condition in stressed_conditionLst:
    myQTp_book[condition]={}
    for Qvalue in [0.001,0.005,0.01,0.02,0.05,0.1]:
        Ttest_pValue=gene_QTp_valuerevDic[condition](Qvalue)
        myQTp_book[condition][Qvalue]=Ttest_pValue
    myQTp_DF=pd.DataFrame(myQTp_book)
    myQTp_DF['Qvalue']=myQTp_DF.index
    myQTp_DF.to_csv('%s_results/%s_Pvalue_Qvalue/Qvalue_TtestPvalue.txt'%(prefix,prefix),index=None,sep='\t')

# store the T test P value cutoff for several commonly used Q value
myQFPR_book={}
for condition in stressed_conditionLst:
    myQFPR_book[condition]={}
    for Qvalue in [0.001,0.005,0.01,0.02,0.05,0.1]:
        FPRvalue=gene_QFPR_valuerevDic[condition](Qvalue)
        myQFPR_book[condition][Qvalue]=FPRvalue
    myQFPR_DF=pd.DataFrame(myQFPR_book)
    myQFPR_DF['Qvalue']=myQFPR_DF.index
    myQFPR_DF.to_csv('%s_results/%s_Pvalue_Qvalue/Qvalue_scoreFPR.txt'%(prefix,prefix),index=None,sep='\t')

# ==========================================
# Calculate Q value for each gene
for condition in stressed_conditionLst:
    QvalueLst=[]; FDRvalueLst=[]
    for gene in processed_geneLst:
        if gene!='0':
            sgRNAnumber=Result_Dic[condition][gene]['sgRNAnumber']
            if sgRNAnumber!=0:
                Ttest_pValue=Result_Dic[condition][gene]['PvalueTtest']
                score=Result_Dic[condition][gene]['score']
                FPRvalue=Result_Dic[condition][gene]['FPRvalue']
                Qvalue=gene_QTp_valueDic[condition](Ttest_pValue)
                FDRvalue=gene_QFPR_valueDic[condition](FPRvalue)
                Result_Dic[condition][gene]['QTpvalue']=Qvalue
                Result_Dic[condition][gene]['FDRvalue']=FDRvalue
                QvalueLst.append(Qvalue); FDRvalueLst.append(FDRvalue)
            else:
                Result_Dic[condition][gene]['QTpvalue']=1
                Result_Dic[condition][gene]['FDRvalue']=1
pickle.dump(Result_Dic,open('%s_results/%s_Result_geneDic.pickle'%(prefix,prefix),'wb'))  #

# ///////////////////////////////////////////////////////////////////
# output all results into a csv file
os.system('mkdir -p %s_results/%s_gene_statistics/'%(prefix,prefix))
for condition in stressed_conditionLst:
    os.system('cat /dev/null > %s_results/%s_gene_statistics/%s_gene_statistics.txt'%(prefix,prefix,condition))
    g=open('%s_results/%s_gene_statistics/%s_gene_statistics.txt'%(prefix,prefix,condition),'r+')
    g.write('gene\tsgRNAnumber\tMedianRAC\tMedianZ\t-Log10Pvalue_MWUtest\tFDRvalue\tFPRvalue\t-Log10Pvalue_Ttest\tQvalue_Ttest\n')
    for gene in processed_geneLst:
        if gene!='0':
            medianRAC=Result_Dic[condition][gene]['MedianRAC']
            medianZ=Result_Dic[condition][gene]['MedianZ']
            sgRNAnumber=Result_Dic[condition][gene]['sgRNAnumber']
            Log10P_MWUtest=Result_Dic[condition][gene]['Log10PvalueMWU']
            FDRvalue=Result_Dic[condition][gene]['FDRvalue']
            FPRvalue=Result_Dic[condition][gene]['FPRvalue']
            Ttest_pValue=Result_Dic[condition][gene]['PvalueTtest']
            Log10Pvalue_Ttest=math.log(Ttest_pValue,10)
            Qvalue=Result_Dic[condition][gene]['QTpvalue']
            g.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(gene,str(sgRNAnumber),str(medianRAC),str(medianZ),str(-1*Log10P_MWUtest),str(FDRvalue),str(FPRvalue),str(-1*Log10Pvalue_Ttest),str(Qvalue)))

    g.close()
