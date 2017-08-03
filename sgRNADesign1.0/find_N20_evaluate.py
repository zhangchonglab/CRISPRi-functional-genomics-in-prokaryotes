# this script is sued to find the N20 site for Cas9 recognition.
# principle:
# find a N20 just upstream of NGG or NAG, 
# Meanwhile, a search against genome is conducted by standalone seqmap tool to remove those with potential offtarget effect

# input: the target nulc seq in fasta format
# input: genome nucl seq in fasta format

# whether this target is within the genome or not
# input: endogenous or exogenous

import os
import sys
import random
from random import randint
from random import choice
import numpy as np

gene_name=sys.argv[1]
gene_seq=sys.argv[2].upper()
genome=sys.argv[3]
N20_number=int(sys.argv[4])
ORFcutoff=float(sys.argv[5]) # within which to design as much as possible
GCcontentMin=int(sys.argv[6])
GCcontentMax=int(sys.argv[7])
strand=str(sys.argv[8])
off_threshold=int(sys.argv[9])
spacer_length=int(sys.argv[10])
PAM=sys.argv[11]
print(strand)
# obtain the reverse complement sequence
def reverse_comp_get(plus_seq):
    CODE={'A':'T','T':'A','C':'G','G':'C','N':'N'} 
    minus_seq=''
    for c in plus_seq:
        minus_seq=minus_seq+CODE[c]
        reverse_comp_seq=minus_seq[::-1]
    return reverse_comp_seq

# obtain the GC content of N20
def GCcontent_extract(seq):
    seqLen=len(seq)
    GCnumber=0
    for nulc in seq:
        if nulc=='G' or nulc=='C':
            GCnumber+=1
    GContentPerc=100*float(GCnumber)/float(seqLen)
    return GContentPerc

# to identify the position of mismatch of two oligos,N20NGG
# N of NGG is defined as zero position GG is defined as -1 and -2, respectively
# the input should be in a N20NGG format, the reverse complement format is strictly forbidden
def mismatch_position_identify(seqOn,seqOff):
    mismatchLst=[]
    seqOnRC=reverse_comp_get(seqOn)
    seqOffRC=reverse_comp_get(seqOff)
    lenOn=len(seqOnRC)
    lenOff=len(seqOffRC)
    if lenOn!=lenOff:
        return 'Incorrect Input of mismatch comparison!'
    for i in range(lenOn):
        if seqOnRC[i]!=seqOffRC[i]:
            mismatchLst.append(i-2)
    return mismatchLst
# thus, the output of this function is a list of position for mismatch with the terminal closed to that of PAM site as 1 position 
# mismatch within NGG is thus -2, -1 and 0, will not be considered in the following function

# scoring system to test the specificity of N20
# according to the penalty and Zhang Feng paper, 2013 NBT, 31 827-32, threshold as <11
# the input should be in a N20NGG format, the reverse complement format is strictly forbidden
def N20_specificity_penalty(seqOn,seqOff):
    penalty=0
    if seqOff[-2:]=='GG':
        mismatchLst=mismatch_position_identify(seqOn,seqOff)
        for mismatch in mismatchLst:
            if mismatch>0 and mismatch<=7:
                penalty+=8
            elif mismatch>7 and mismatch<=12:
                penalty+=4.5
            elif mismatch>12:
                penalty+=2.5
    elif seqOff[-2:]=='AG':
        mismatchLst=mismatch_position_identify(seqOn,seqOff)
        for mismatch in mismatchLst:
            if mismatch>0 and mismatch<=7:
                penalty+=10
            elif mismatch>7 and mismatch<=12:
                penalty+=7
            elif mismatch>12:
                penalty+=3
    else:
        penalty=40
    return penalty

# the function to test whether one sgRNA has off-target effect above threshold or not
# the input is the list of penalty [0,0,...]
def sgRNA_judge(offtarget_penalty,off_threshold):
    offtarget_penalty.sort()
    if len(offtarget_penalty)<2:
        return False
    elif len(offtarget_penalty)==2 and offtarget_penalty[0]==0 and offtarget_penalty[1]==0:
        return True

    elif len(offtarget_penalty)>2 and offtarget_penalty[2]>off_threshold: # the best offtarget has big penalty
        return True # this sgRNA is ok
    else:
        return False

# define a function to calculate the average and stdev of spacing for a sorted list
def mean_stdev_lst(ID_lst):
    ID_lst.sort()
    spacing_lst=[]
    for i,item in enumerate(ID_lst):
        if i!=0:
            spacing_lst.append(float(ID_lst[i]-ID_lst[i-1]))
    return np.mean(spacing_lst),np.std(spacing_lst)

# define a function to extract ID evenly distriuted along CDS
def even_distribute(ID_lst,wanted_number):
    ID_lst.sort()
    if wanted_number>=len(ID_lst):
        return ID_lst
    elif wanted_number<len(ID_lst) and wanted_number==1:
        random.shuffle(ID_lst)
        returnID_lst=ID_lst[:wanted_number]
        return returnID_lst
    else:
        perfect_spacing=float(ID_lst[-1]-ID_lst[0])/float(wanted_number-1)
        random.shuffle(ID_lst)
        returnID_lst=ID_lst[:wanted_number]
        spacing_mean=0
        spacing_stdev=10000
        spacing_mean_record=0
        spacing_stdev_record=10000
        for counter in range(1000):
            random_subset=random.sample(range(0,len(ID_lst)),wanted_number)
            test_subset=[]
            for i,ID in enumerate(ID_lst):
                if i in random_subset:
                    test_subset.append(ID)
            spacing_mean,spacing_stdev=mean_stdev_lst(test_subset)
            if spacing_mean>0.8*spacing_mean_record and spacing_stdev<spacing_stdev_record:
                returnID_lst=test_subset
                spacing_stdev_record=spacing_stdev
                spacing_mean_record=spacing_mean
        if spacing_mean_record>0.9*perfect_spacing and spacing_stdev_record<0.2*perfect_spacing:
            return returnID_lst
        else:
            for counter in range(100000):
                random_subset=random.sample(range(0,len(ID_lst)),wanted_number)
                test_subset=[]
                for i,ID in enumerate(ID_lst):
                    if i in random_subset:
                        test_subset.append(ID)
                spacing_mean,spacing_stdev=mean_stdev_lst(test_subset)
                if spacing_mean>0.95*spacing_mean_record and spacing_stdev<spacing_stdev_record:
                    returnID_lst=test_subset
                    spacing_stdev_record=spacing_stdev
                    spacing_mean_record=spacing_mean
                if spacing_mean_record>0.9*perfect_spacing and spacing_stdev_record<0.2*perfect_spacing:
                    return returnID_lst
            return returnID_lst
# /////////////////////////////////////////////////////////////////////////////////////////////
# dic to store all the candidate N20NGG
# {'ID':seq, ...}
Candidate_N20NGGDic={}

# {'ID':total_number, ...} including target, thus 1 is ok for selection
Candidate_N20mismatchDic={}


# put all the potential N20NGG into a fasta file and a dictionary with N20NGG as index
# candidateoff.fasta to store all the N20-N(abbr)G(A)G seq in fasta format
os.system('cat /dev/null > N20%sofftarget.fasta'%(PAM))
sgRNA_ID_collection=[]
sgRNA_existence_flag=False
for position in range(len(gene_seq)):
    geneLen=len(gene_seq)
    if strand=='nontemplate':
        if position<(geneLen-spacer_length-len(PAM)+1):
            if gene_seq[position:(position+len(PAM)-1)]=='CC':
                sgRNA_existence_flag=True
                N20=reverse_comp_get(gene_seq[(position+len(PAM)):(position+spacer_length+len(PAM))])
                candidateN20NGG=reverse_comp_get(gene_seq[position:(position+spacer_length+len(PAM))])
                sgRNA_ID=gene_name+'_'+str(position+1)
                if sgRNA_ID not in sgRNA_ID_collection:
                    sgRNA_ID_collection.append(int(sgRNA_ID.split('_')[1]))
                Candidate_N20NGGDic[sgRNA_ID]=candidateN20NGG
                Candidate_N20mismatchDic[sgRNA_ID]=[]
                flag=True  # monitor whether no off-target site idnetified
                target1Nameline='>%s-Target1'%(sgRNA_ID)
                potentialOffTarget1=N20+'NGG'
                target2Nameline='>%s-Target2'%(sgRNA_ID)
                potentialOffTarget2=N20+'NAG'
                os.system('echo "%s\n%s\n%s\n%s\n" >> N20%sofftarget.fasta'%(target1Nameline,potentialOffTarget1,target2Nameline,potentialOffTarget2,PAM))
    elif strand=='template':
        if position>=(spacer_length+len(PAM)-1):
            if gene_seq[(position-1):position+1]=='GG':
                sgRNA_existence_flag=True
                N20=gene_seq[(position-spacer_length-len(PAM)+1):(position-len(PAM)+1)]
                candidateN20NGG=gene_seq[(position-spacer_length-len(PAM)+1):position+1]
                sgRNA_ID=gene_name+'_'+str(position-1-spacer_length)
                if sgRNA_ID not in sgRNA_ID_collection:
                    sgRNA_ID_collection.append(int(sgRNA_ID.split('_')[1]))
                Candidate_N20NGGDic[sgRNA_ID]=candidateN20NGG
                Candidate_N20mismatchDic[sgRNA_ID]=[]
                flag=True  # monitor whether no off-target site idnetified
                target1Nameline='>%s-Target1'%(sgRNA_ID)
                potentialOffTarget1=N20+'NGG'
                target2Nameline='>%s-Target2'%(sgRNA_ID)
                potentialOffTarget2=N20+'NAG'
                os.system('echo "%s\n%s\n%s\n%s\n" >> N20%sofftarget.fasta'%(target1Nameline,potentialOffTarget1,target2Nameline,potentialOffTarget2,PAM))
# After this step, all the potential N20NG(A)G are involved into N20NGGofftarget.fasta and subjected to seqmap for analysis

# //////////////////////
# write in sgRNAs into the file with user defined method
# candidateN20.fasta to store all the N20-N(exact)G(A)G seq in fasta format
os.system('cat /dev/null > %s_candidateN20.fasta'%(gene_name))
os.system('cat /dev/null > %s_candidateN20NGG.fasta'%(gene_name))
os.system('cat /dev/null > %s_statistics.txt'%(gene_name))
os.system('cat /dev/null > %s_sgRNA_statistics.txt'%(gene_name))
g=open('%s_candidateN20.fasta'%(gene_name),'r+')
k=open('%s_candidateN20NGG.fasta'%(gene_name),'r+')
m=open('%s_statistics.txt'%(gene_name),'r+')
n=open('%s_sgRNA_statistics.txt'%(gene_name),'r+')
# there is at least one sgRNA
if sgRNA_existence_flag:
    # conduct seqmap to find the potential offtarget
    os.system('./seqmap 5 N20%sofftarget.fasta %s seqmapOut.temp.dat /output_all_matches'%(PAM,genome))

    f=open('seqmapOut.temp.dat','r')

    # to find the mismatch hit in genome for each sgRNA
    for line in f:
        row=line.rstrip().split('\t')
        if row[0]!='trans_id':
            sgRNA_ID=row[3].split('-')[0]
            guide_seq=row[4]
            mismatch_seq=row[2]
            penalty=N20_specificity_penalty(guide_seq,mismatch_seq) 
            Candidate_N20mismatchDic[sgRNA_ID].append(penalty)
    f.close()

    os.system('rm seqmapOut*dat N20%sofftarget.fasta'%(PAM))

    # store all sgRNAs passing the quality control
    qualitified_sgRNAID_list=[]
    qualitified_sgRNAID_dic={}
    for ID in sgRNA_ID_collection:
        sgRNA_ID='%s_%d'%(gene_name,ID)
        if sgRNA_judge(Candidate_N20mismatchDic[sgRNA_ID],off_threshold):  # pass the first threshold
            N20NGG=Candidate_N20NGGDic[sgRNA_ID]
            N20=Candidate_N20NGGDic[sgRNA_ID][:spacer_length]
            if GCcontent_extract(N20NGG)<GCcontentMax and GCcontent_extract(N20NGG)>GCcontentMin:  # GC content ok
                qualitified_sgRNAID_dic[sgRNA_ID]=[N20,N20NGG]
                qualitified_sgRNAID_list.append(ID)
    # sort the list to get the 5 proximal firstly
    qualitified_sgRNAID_list.sort()

    # record sgRNA number:
    number_record=0
    # to put the sgRNA passing the threshold into the doc
    # sgRNA within the defined region of CDS
    regionThre=float(len(gene_seq)*ORFcutoff) # within which to design as many as possible
    sgRNAIDoutside=[]
    for ID in qualitified_sgRNAID_list:
        sgRNA_ID='%s_%d'%(gene_name,ID)
        if ID<regionThre:
            sgRNA_ID='%s_%d'%(gene_name,ID)
            N20=qualitified_sgRNAID_dic[sgRNA_ID][0]
            N20NGG=qualitified_sgRNAID_dic[sgRNA_ID][1]
            g.write('>%s\n%s\n'%(sgRNA_ID,N20))
            k.write('>%s\n%s\n'%(sgRNA_ID,N20NGG))
            n.write('%s\t%s\t%s\n'%(sgRNA_ID,str(float(ID)/float(len(gene_seq))),str(GCcontent_extract(N20NGG))))
            number_record+=1
        else:
            sgRNAIDoutside.append(ID)
        if number_record>=N20_number:
            break

    # to boost designed sgRNA to N20_number
    if number_record<N20_number:
        wanted_number=N20_number-number_record
        evenly_distributedID_lst=even_distribute(sgRNAIDoutside,wanted_number)
        for ID in evenly_distributedID_lst:
            sgRNA_ID='%s_%d'%(gene_name,ID)
            N20=qualitified_sgRNAID_dic[sgRNA_ID][0]
            N20NGG=qualitified_sgRNAID_dic[sgRNA_ID][1]
            g.write('>%s\n%s\n'%(sgRNA_ID,N20))
            k.write('>%s\n%s\n'%(sgRNA_ID,N20NGG))
            n.write('%s\t%s\t%s\n'%(sgRNA_ID,str(float(ID)/float(len(gene_seq))),str(GCcontent_extract(N20NGG))))
            number_record+=1
            if number_record>=N20_number:
                break
    print '%d sgRNA detected for %s'%(number_record,gene_name)
    m.write('%s\t%d\t%d\n'%(gene_name,len(gene_seq),number_record))

# no sgRNA found for this gene
else:
    print 'no sgRNA detected for %s'%(gene_name)
    m.write('%s\t%d\t0\n'%(gene_name,len(gene_seq)))
print ''
print ''
k.close()
g.close()
m.close()
n.close()
