import os
import sys
import numpy as np
from scipy.stats import norm
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from matplotlib import cm
import pickle
from scipy.stats.mstats import gmean
import math
import warnings

warnings.simplefilter("error")

exp_configure=sys.argv[1]
name_configure=sys.argv[2]
normalizedData=sys.argv[3]
# method to set the control sgRNA set: 'all or NC'
control_setting=sys.argv[4]
if control_setting not in ['all','NC']:
    print 'incorrect control setting!'
    sys.exit()
prefix=sys.argv[5]

# number of reads in initial library that below this threshold won't be incoperated in to the following analysis pipleine 
readsThreshold=float(sys.argv[6])     #need choose 
os.system('mkdir %s_results/'%(prefix))

# ///////////////////////////////
# process the exp_configure data
LibraryLst=[]
ConditionLst=[]
ConditionLibrary={}
stress_control_pair_number=1
Library_readsnumber={}
f=open(exp_configure,'r')
for i,line in enumerate(f):
    row=line.rstrip().split('\t')
    if i==0:
        ConditionLst=row[1:]
        stress_control_pair_number=(len(ConditionLst)-1)/2
        for condition in ConditionLst:
             ConditionLibrary[condition]=[]
    else:
        library=row[0]
        Library_readsnumber[library]=0
        LibraryLst.append(library)
        condition_index=[i for i,j in enumerate(row[1:]) if j=='1']
        conditions=[ConditionLst[i] for i in condition_index]
        for condition in conditions:
            ConditionLibrary[condition].append(library)
f.close()

# stressed condition list that need to be processed into final files
stressed_conditionLst=[]
for pair_number in range(stress_control_pair_number):
    condition='stress'+str(pair_number+1)
    stressed_conditionLst.append(condition)

# stressed condition label and name dictionary
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

stressed_conditionDic={}
stressed_condition_nameLst=[]
f=open(name_configure,'r')
for i,line in enumerate(f):
    row=line.rstrip().split('\t')
    condition=row[0]
    name=row[1]
    if 'stress'==condition[:6] and RepresentsInt(condition[6:]):
        stressed_conditionDic[condition]=name
        stressed_condition_nameLst.append(name)
    else:
        print ('incorrect name configure setting, check line %d in name configure file'%(i+1))
        sys.exit(-1)
f.close()
if len(stressed_condition_nameLst)!=len(stressed_conditionLst):
    print ('inconsistent setting between name and experiment configure files! different number of experiments detected in these two files!')
    sys.exit(-1)
pickle.dump(stressed_condition_nameLst,open('%s_results/%s_stressed_conditionLst.pickle'%(prefix,prefix), 'wb'))

# ////////////////////////////////////////////////////
# calculate the total reads of each library
f=open(normalizedData,'r')
for i,line in enumerate(f):
    row=line.rstrip().split('\t')
    if i==0:
        LibraryKeys=row[2:]
    else:
        for ii,reads in enumerate(row[2:]):
            Library_readsnumber[LibraryKeys[ii]]+=float(reads)
f.close()

# //////////////////////
# process the normalized abundance data
# sgRNA: gene,LibAbundanceDic:{Lib1:xx,Lib2:xx, ..} LibReadsDic:{Lib1:xx,Lib2:xx, ..},ConditionAbundanceDic:{condition1:xx,condition:xx,..},ConditionReadsDic:{condition1:xx,condition:xx,..}

removed_sgRNALst=[]
processed_sgRNALst=[]
Processed_sgRNADic={}
processed_geneLst=[]
f=open(normalizedData,'r')
for i,line in enumerate(f):
    row=line.rstrip().split('\t')
    if i==0:
        LibraryKeys=row[2:]
    else:
        sgRNA=row[0]
        gene=row[1]
        Lib_abundanceDic=dict.fromkeys(LibraryKeys)
        Lib_readsDic=dict.fromkeys(LibraryKeys)
        Lib_abundance_vs_initialDic=dict.fromkeys(LibraryKeys)
        Condition_abundanceDic=dict.fromkeys(ConditionLst)
        Condition_readsDic=dict.fromkeys(ConditionLst)
        Condition_relative_deviationDic=dict.fromkeys(ConditionLst)
        Condition_abundance_vs_initialDic=dict.fromkeys(ConditionLst)
        for ii,reads in enumerate(row[2:]):
            if float(reads)<=1:
                reads=1
            Lib_abundanceDic[LibraryKeys[ii]]=math.log(float(reads)/Library_readsnumber[LibraryKeys[ii]],2)
            Lib_readsDic[LibraryKeys[ii]]=float(reads)
        # combine the biological replicates by geometric average
        for condition in ConditionLst:
            librariesOfThisCondition=ConditionLibrary[condition]
            abundanceOfThisCondition=[]
            readsOfThisCondition=[]
            for library in librariesOfThisCondition:
                abundanceOfThisCondition.append(Lib_abundanceDic[library])
                readsOfThisCondition.append(Lib_readsDic[library])
            # abundance is the log2 value, thus log2 genometric average is the arithmetic average
            conditionAbundance=np.mean(abundanceOfThisCondition)
            conditionReads=gmean(readsOfThisCondition)
            conditionRelativeDeviation=np.std(readsOfThisCondition)/np.mean(readsOfThisCondition)
            Condition_abundanceDic[condition]=conditionAbundance
            Condition_readsDic[condition]=conditionReads
            Condition_relative_deviationDic[condition]=conditionRelativeDeviation
        # judge whether to subject to analysis or remove
        if Condition_readsDic['initial']>=readsThreshold:
            processed_sgRNALst.append(sgRNA)
            if gene not in processed_geneLst:
                processed_geneLst.append(gene)
            Condition_abundance_vs_initial=Condition_abundanceDic['initial']
            Lib_abundance_vs_initial=Condition_abundanceDic['initial']
            for library in LibraryLst:
                Lib_abundance_vs_initialDic[library]=Lib_abundanceDic[library]-Lib_abundance_vs_initial
            for condition in ConditionLst:
                Condition_abundance_vs_initialDic[condition]=Condition_abundanceDic[condition]-Condition_abundance_vs_initial
            Processed_sgRNADic[sgRNA]={}
            Processed_sgRNADic[sgRNA]['Gene']=gene
            Processed_sgRNADic[sgRNA]['LibAbundanceDic']=Lib_abundanceDic
            Processed_sgRNADic[sgRNA]['LibReadsDic']=Lib_readsDic
            Processed_sgRNADic[sgRNA]['LibAbundance_vs_initialDic']=Lib_abundance_vs_initialDic
            Processed_sgRNADic[sgRNA]['ConditionAbundanceDic']=Condition_abundanceDic
            Processed_sgRNADic[sgRNA]['ConditionReadsDic']=Condition_readsDic
            Processed_sgRNADic[sgRNA]['ConditionAbundance_vs_initialDic']=Condition_abundance_vs_initialDic
            Processed_sgRNADic[sgRNA]['ConditionRelative_deviationDic']=Condition_relative_deviationDic
        else:
            removed_sgRNALst.append(sgRNA)
f.close()

# calculate the relative abundance change between correspondingly stressed and control conditions
for sgRNA in processed_sgRNALst:
    Processed_sgRNADic[sgRNA]['RelativeAbundanceChange']={}
    for condition in stressed_conditionLst:
        pair_number=condition.split('stress')[1]
        abundance_stress=Processed_sgRNADic[sgRNA]['ConditionAbundanceDic']['stress%s'%(pair_number)]
        abundance_control=Processed_sgRNADic[sgRNA]['ConditionAbundanceDic']['control%s'%(pair_number)]
        Processed_sgRNADic[sgRNA]['RelativeAbundanceChange'][condition]=abundance_stress-abundance_control



# /////////////////////////////
# write removed sgRNAs
os.system('cat /dev/null > %s_results/%s.removed.sgRNA.txt'%(prefix,prefix))
g=open('%s_results/%s.removed.sgRNA.txt'%(prefix,prefix),'r+')
for sgRNA in removed_sgRNALst:
    g.write(sgRNA+'\n')
g.close()
os.system('mkdir %s_results/removed.sgRNA/'%(prefix))
os.system('mv %s_results/%s.removed.sgRNA.txt %s_results/removed.sgRNA/'%(prefix,prefix,prefix))

# ////////////////////////////////////////////
# give replicates consistence data
os.system('cat /dev/null > %s_results/%s_replicates_reads_statistics.txt'%(prefix,prefix))
k=open('%s_results/%s_replicates_reads_statistics.txt'%(prefix,prefix),'r+')
k.write('condition\tpearson correlation coefficient\tP value\n')
for condition in ConditionLst:
    if len(ConditionLibrary[condition])>=2: # replicates existence check
        os.system('cat /dev/null > %s_results/%s_%s_replicates.txt'%(prefix,prefix,condition))
        g=open('%s_results/%s_%s_replicates.txt'%(prefix,prefix,condition),'r+')
        writtenLine='sgRNA\t'
        for library in ConditionLibrary[condition]:
            writtenLine+='%s_abundance\t%s_reads\t%s_abundance_vs_initial\t'%(library,library,library)
        g.write(writtenLine[:-1]+'\n')
        Xaxis=[]
        Yaxis=[]
        for sgRNA in processed_sgRNALst:
            record=1
            writtenLine=''
            for library in ConditionLibrary[condition]:
                abundance=Processed_sgRNADic[sgRNA]['LibAbundanceDic'][library]
                reads=Processed_sgRNADic[sgRNA]['LibReadsDic'][library]
                abundance_vs_initial=Processed_sgRNADic[sgRNA]['LibAbundance_vs_initialDic'][library]
                if record==1:
                    writtenLine+='%s\t'%(sgRNA)
                    Xlibrary=library
                    Xaxis.append(reads)
                elif record==2:
                    Yaxis.append(reads)
                    Ylibrary=library
                writtenLine+='%s\t%s\t%s\t'%(abundance,reads,abundance_vs_initial)
                record+=1
            g.write(writtenLine[:-1]+'\n')
        pearsonCE,pValue=pearsonr(np.array(Xaxis),np.array(Yaxis))
        xy = np.vstack([np.array(Xaxis),np.array(Yaxis)])
        z = gaussian_kde(xy)(xy)
        plt.scatter(np.array(Xaxis),np.array(Yaxis),c=z, s=10, edgecolor='none')
        plt.xlabel('%s normalized reads'%(Xlibrary),fontsize=16)
        plt.ylabel('%s normalized reads'%(Ylibrary),fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim((0.1,max(Xaxis)*3))
        plt.ylim((0.1,max(Yaxis)*3))
        plt.text(0.15,max(Yaxis),'Pearson Correlation Coefficient=%s'%(("{0:.3f}".format(pearsonCE))),fontsize=14)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.95)
        plt.savefig('%s_results/%s_%s_replicates.png'%(prefix,prefix,condition),dpi=400)
        plt.clf()
        k.write('%s\t%s\t%s\n'%(condition,pearsonCE,pValue))
        g.close()
k.close()

# move these flat files into a sub directory named replicate_consistence
os.system('mkdir %s_results/replicate_consistence/'%(prefix))
os.system('mv %s_results/*_replicates*txt %s_results/replicate_consistence/'%(prefix,prefix))
os.system('mv %s_results/*_replicates.png %s_results/replicate_consistence/'%(prefix,prefix))


# /////////////////////////////////////////////////
# negative control sgRNA normal distribution fitting to ND
NCsgRNA_dic={}
NCsgRNA_ND_dic={}
os.system('cat /dev/null > %s_results/%s_NCsgRNA_ND.txt'%(prefix,prefix))
g=open('%s_results/%s_NCsgRNA_ND.txt'%(prefix,prefix),'r+')
g.write('condition\tmedian\tmean\tstdev\n')
for condition in stressed_conditionLst:
    condition_name=stressed_conditionDic[condition]
    NCsgRNA_dic[condition]=[]
    NCsgRNA_ND_dic[condition]={}
    for sgRNA in processed_sgRNALst:
        if control_setting=='NC': # use negative control sgRNA set for normalization
            if Processed_sgRNADic[sgRNA]['Gene']=='0':
                relative_abundance=Processed_sgRNADic[sgRNA]['RelativeAbundanceChange'][condition]
                NCsgRNA_dic[condition].append(relative_abundance)
        elif control_setting=='all':
            relative_abundance=Processed_sgRNADic[sgRNA]['RelativeAbundanceChange'][condition]
            NCsgRNA_dic[condition].append(relative_abundance)
    meanZ, stdZ = norm.fit(NCsgRNA_dic[condition])
    medianZ=np.median(np.array(NCsgRNA_dic[condition]))
    g.write('%s\t%s\t%s\t%s\n'%(condition_name,str(medianZ),str(meanZ),str(stdZ)))
    NCsgRNA_ND_dic[condition]['Mean']=meanZ
    NCsgRNA_ND_dic[condition]['Median']=medianZ
    NCsgRNA_ND_dic[condition]['Stdev']=stdZ
g.close()

# //////////////////////////////////////////////////
# start to normalize sgRNA relative abundance by NC sgRNA ND medianZ and stdZ
for sgRNA in processed_sgRNALst:
    Processed_sgRNADic[sgRNA]['Normalized_RelativeAbundanceChange']={}
    Processed_sgRNADic[sgRNA]['Zscore']={}
    for condition in stressed_conditionLst:
        relativeAbundanceChange=Processed_sgRNADic[sgRNA]['RelativeAbundanceChange'][condition]
        medianZ_NC=NCsgRNA_ND_dic[condition]['Median']
        stdZ=NCsgRNA_ND_dic[condition]['Stdev']
        normalized_change=(relativeAbundanceChange-medianZ_NC)
        Zscore=(relativeAbundanceChange-medianZ_NC)/stdZ
        Processed_sgRNADic[sgRNA]['Normalized_RelativeAbundanceChange'][condition]=normalized_change
        Processed_sgRNADic[sgRNA]['Zscore'][condition]=Zscore

# ////////////////////////////////////////////////////////
# file for further processing in gene level calculation, store relative abundance change
# {gene1:{stressed_cpn1:{sgRNA1:xx, sgRNA2:xx, ..},stressed_con2:[], ..}, gene2:{}, ..}
# Here xx refers to the relative abundance change between stressed and corresponding control conditions
Processed_geneDic={gene:{condition_name:{} for condition_name in stressed_condition_nameLst} for gene in processed_geneLst}
# file for further processing in gene level calculation, store Zscore
# {gene1:{stressed_cpn1:{sgRNA1:xx, sgRNA2:xx, ..},stressed_con2:[], ..}, gene2:{}, ..}
# Here xx refers to the Zscore obtained by relative abundance change/stdZ between stressed and corresponding control conditions
Processed_geneZDic={gene:{condition_name:{} for condition_name in stressed_condition_nameLst} for gene in processed_geneLst}

for sgRNA in processed_sgRNALst:
    gene=Processed_sgRNADic[sgRNA]['Gene']
    for condition in stressed_conditionLst:
        condition_name=stressed_conditionDic[condition]
        index=condition.split('stress')[1]
        relevant_control_condition='control%s'%(index)
        control_read=Processed_sgRNADic[sgRNA]['ConditionReadsDic'][relevant_control_condition]
        if control_read>=readsThreshold:
            relative_abundance_change=Processed_sgRNADic[sgRNA]['Normalized_RelativeAbundanceChange'][condition]
            Processed_geneDic[gene][condition_name][sgRNA]=relative_abundance_change
            Zscore=Processed_sgRNADic[sgRNA]['Zscore'][condition]
            Processed_geneZDic[gene][condition_name][sgRNA]=Zscore

pickle.dump(Processed_geneDic,open('%s_results/%s_Processed_geneDic.pickle'%(prefix,prefix), 'wb'))
pickle.dump(Processed_geneZDic,open('%s_results/%s_Processed_geneZDic.pickle'%(prefix,prefix), 'wb'))
pickle.dump(processed_geneLst,open('%s_results/%s_processed_geneLst.pickle'%(prefix,prefix), 'wb'))
pickle.dump(processed_sgRNALst,open('%s_results/%s_processed_sgRNALst.pickle'%(prefix,prefix), 'wb'))

# //////////////////////////////////////////////
# draw the NC sgRNA distribution in each condition
os.system('cat /dev/null > %s_results/%s_NCsgRNA_normalized_ND.txt'%(prefix,prefix))
g=open('%s_results/%s_NCsgRNA_normalized_ND.txt'%(prefix,prefix),'r+')
g.write('condition\tmedian\tmean\tstdev\n')
NCsgRNA_dic={}
cm_subsection = np.linspace(0.25, 0.35, len(stressed_conditionLst))
colors = [ cm.jet(x) for x in cm_subsection ]
for i,condition in enumerate(stressed_conditionLst):
    condition_name=stressed_conditionDic[condition]
    NCsgRNA_dic[condition]=[]
    for sgRNA in processed_sgRNALst:
        if control_setting=='NC': # use negative control sgRNA set for normalization
            if Processed_sgRNADic[sgRNA]['Gene']=='0':
                relative_abundance=Processed_sgRNADic[sgRNA]['Normalized_RelativeAbundanceChange'][condition]
                NCsgRNA_dic[condition].append(relative_abundance)
        elif control_setting=='all':
            relative_abundance=Processed_sgRNADic[sgRNA]['Normalized_RelativeAbundanceChange'][condition]
            NCsgRNA_dic[condition].append(relative_abundance)
    meanZ, stdZ = norm.fit(NCsgRNA_dic[condition])
    data=np.array(NCsgRNA_dic[condition])
    bins=np.linspace(np.nanmin(data),np.nanmax(data),20)
    plt.hist(data,bins,color=colors[i])
    plt.xlabel('Log2 (normalzied sgRNA relative abundance change)',fontsize=14)
    plt.ylabel('sgRNA number',fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.95)
    plt.savefig('%s_results/%s_%s_NCsgRNAND.png'%(prefix,prefix,condition_name),dpi=400)
    plt.clf()
    medianZ=np.median(np.array(NCsgRNA_dic[condition]))
    g.write('%s\t%s\t%s\t%s\n'%(condition_name,str(medianZ),str(meanZ),str(stdZ)))
g.close()
os.system('mkdir %s_results/NCsgRNA_ND/'%(prefix))
os.system('mv %s_results/%s_NCsgRNA*.txt %s_results/*.png %s_results/NCsgRNA_ND/'%(prefix,prefix,prefix,prefix))

# /////////////////////////////////////////////////////
# start to output all sgRNA information
os.system('mkdir %s_results/%s_sgRNA_statistics/'%(prefix,prefix))
os.system('mkdir %s_results/%s_sgRNA_statistics/library_level/'%(prefix,prefix))
for library in LibraryLst:
    os.system('cat /dev/null > %s_results/%s_sgRNA_statistics/library_level/%s_sgRNA_%s_results.txt'%(prefix,prefix,prefix,library))
    g=open('%s_results/%s_sgRNA_statistics/library_level/%s_sgRNA_%s_results.txt'%(prefix,prefix,prefix,library),'r+')
    g.write('sgRNA\tgene\t%s_Log2_abundnace\t%s_reads\t%s_Log2_abundnace_vs_initial\n'%(library,library,library))
    for sgRNA in Processed_sgRNADic:
        gene=Processed_sgRNADic[sgRNA]['Gene']
        abundance=Processed_sgRNADic[sgRNA]['LibAbundanceDic'][library]
        reads=Processed_sgRNADic[sgRNA]['LibReadsDic'][library]
        abundance_vs_initial=Processed_sgRNADic[sgRNA]['LibAbundance_vs_initialDic'][library]
        g.write('%s\t%s\t%s\t%s\t%s\n'%(sgRNA,gene,str(abundance),str(reads),str(abundance_vs_initial)))
    g.close()

os.system('mkdir %s_results/%s_sgRNA_statistics/condition_level/'%(prefix,prefix))
for condition in ConditionLst:
    os.system('cat /dev/null > %s_results/%s_sgRNA_statistics/condition_level/%s_sgRNA_%s_results.txt'%(prefix,prefix,prefix,condition))
    g=open('%s_results/%s_sgRNA_statistics/condition_level/%s_sgRNA_%s_results.txt'%(prefix,prefix,prefix,condition),'r+')
    g.write('sgRNA\tgene\t%s_Log2_abundnace\t%s_reads\t%s_Log2_abundnace_vs_initial\t%s_relative_deviation\n'%(condition,condition,condition,condition))
    for sgRNA in Processed_sgRNADic:
        gene=Processed_sgRNADic[sgRNA]['Gene']
        abundance=Processed_sgRNADic[sgRNA]['ConditionAbundanceDic'][condition]
        reads=Processed_sgRNADic[sgRNA]['ConditionReadsDic'][condition]
        abundance_vs_initial=Processed_sgRNADic[sgRNA]['ConditionAbundance_vs_initialDic'][condition]
        relative_deviation=Processed_sgRNADic[sgRNA]['ConditionRelative_deviationDic'][condition]
        g.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(sgRNA,gene,str(abundance),str(reads),str(abundance_vs_initial),str(relative_deviation)))
    g.close()

os.system('mkdir %s_results/%s_sgRNA_statistics/combined_condition_level/'%(prefix,prefix))
for condition in stressed_conditionLst:
    condition_name=stressed_conditionDic[condition]
    os.system('cat /dev/null > %s_results/%s_sgRNA_statistics/combined_condition_level/sgRNA_combined_%s_results.txt'%(prefix,prefix,condition_name))
    g=open('%s_results/%s_sgRNA_statistics/combined_condition_level/sgRNA_combined_%s_results.txt'%(prefix,prefix,condition_name),'r+')
    g.write('sgRNA\tgene\t%s_relative_abundnace_change\t%s_normalized_change\t%s_Zscore\tQuality\n'%(condition_name,condition_name,condition_name))
    index=condition.split('stress')[1]
    relevant_control_condition='control%s'%(index)
    for sgRNA in Processed_sgRNADic:
        gene=Processed_sgRNADic[sgRNA]['Gene']
        relative_abundance_change=Processed_sgRNADic[sgRNA]['RelativeAbundanceChange'][condition]
        normalized_change=Processed_sgRNADic[sgRNA]['Normalized_RelativeAbundanceChange'][condition]
        Zscore=Processed_sgRNADic[sgRNA]['Zscore'][condition]
        control_read=Processed_sgRNADic[sgRNA]['ConditionReadsDic'][relevant_control_condition]
        if control_read>=readsThreshold:
            flag='Good'
        else:
            flag='Bad'
        g.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(sgRNA,gene,str(relative_abundance_change),str(normalized_change),str(Zscore),flag))
    g.close()
