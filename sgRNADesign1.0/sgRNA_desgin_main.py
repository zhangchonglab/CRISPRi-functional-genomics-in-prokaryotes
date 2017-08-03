# this script is used to design N20 for target nucl seq in a high-throughput fashion
# in this case, all the target seqs are the ORF part and thus non-template strand is targeted

# two input files:

# target seq 
# fasta format

# target genome
# fasta format

# N20_number
# a maximum number of N20 for each target gene

# in the following part, for each target nulc seq, son function find_N20_evaluate.py is cited to compute N20 and eliminate those with potential drawbacks (off-target due to suffering mismatch, poly U/T, GC content, sgRNA structure disruption)

# all the ok N20 is saved and returned back to this main function and written into the output file
import os
import sys
import matplotlib.pyplot as plt
import numpy as np 

########################
#process the configure file
########################
configure=open(sys.argv[1],'r') 
itemDic={}
for rec in configure:
    rec=rec.strip().split('\t')
    item=rec[0]
    variable=rec[1]
    itemDic[item]=variable
	

##########################
sgRNA_number=int(itemDic['sgRNA_number'])
####
spacer_length=20
####
if int(itemDic['GCcontentMin'])>100 or int(itemDic['GCcontentMin'])<0:
    print('please input the positive integer for GCcontentMin')
    sys.exit()
GCcontentMin=int(itemDic['GCcontentMin'])
####
if int(itemDic['GCcontentMax'])>100 or int(itemDic['GCcontentMax'])<0:
    print('please input the positive integer for GCcontentMax')
    sys.exit()
GCcontentMax=int(itemDic['GCcontentMax'])
####
if int(itemDic['off_threshold'])<10:
    print('please input the positive integer for off_threshold')
    sys.exit()
off_threshold=int(itemDic['off_threshold'])
####
if float(itemDic['ORFcutoff'])>1 or float(itemDic['ORFcutoff'])<0:
    print('please input the correct decimal for ORFcutoff')
    sys.exit()
ORFcutoff=float(itemDic['ORFcutoff'])
####
PAM='NGG'
####
strand=itemDic['strand']
if strand!='template' and strand!='nontemplate':
    print('please input the correct strand')
    sys.exit()
####
negative=itemDic['negative']
if negative!='yes' and negative!='no':
    print('please input the correct negative')
    sys.exit()
####
negative_number=itemDic['negative_number']
####
indexFile=itemDic['indexFile']   # ptt or rnt files
targetFasta=itemDic['targetFasta']  #ffn or frn files
genome=itemDic['genome'] #fna
####
if itemDic['prefix']=='':
    prefix='output'
else:
    prefix=itemDic['prefix']
##########################

# ///////////////////////////////////////////////////
# processing the ptt or rnt file
genePositionLst=[]
genePositionDic={}
f=open(indexFile,'r')
for i,line in enumerate(f):
    if i>=3:
        row=line.rstrip().split('\t')
        chain=row[1]
        gene=row[4]
        if '_' not in row[5]:
            synonym=row[5]
        else:
            synonym=row[5].split('_')[0]+row[5].split('_')[1]
        if chain=='+':
            start=row[0].split('..')[0]
        if chain=='-':
            start=row[0].split('..')[1]
        genePositionDic[start]=synonym
        genePositionLst.append(start)
f.close()

# to record gene name difference in fasta and index file
os.system('cat /dev/null > %s.gene_name.notconsistent.list.txt'%(prefix))
os.system('cat /dev/null > %s.gene_start.notconsistent.list.txt'%(prefix))
g=open('%s.gene_name.notconsistent.list.txt'%(prefix),'r+')
g.write('index\tfasta\tconsistence\n')
k=open('%s.gene_start.notconsistent.list.txt'%(prefix),'r+')
k.write('fastaGene\tfastaStart\n')
os.system('cat /dev/null > %s.RefSeq_geneID.txt'%(prefix))
m=open('%s.RefSeq_geneID.txt'%(prefix),'r+')
m.write('#RefSeqID\tgeneID\n')
# start to check each entry of target fasta file
f=open(targetFasta,'r')
target_fasta_dic={}
geneID='' # to store the protein ID now
for line in f:
    if line[0]=='>':
        flag=False
        row=line.rstrip().split(' [')
        
        location_item=line.rstrip().split(':')[1].split(' ')[0]
        if 'c' in location_item:
            start=location_item[:-1].split('c')[1].split('-')[0]
           
            if start not in genePositionLst:
                print ('%s,Yes'%(start))
        else:
            start=location_item[:-1].split('-')[0]
        if start in genePositionLst:
            geneID=genePositionDic[start]
            target_fasta_dic[geneID]=''
            flag=True
    elif flag:
        target_fasta_dic[geneID]+=line.rstrip()



f.close()
g.close()
k.close()
m.close()
os.system('rm -r %s'%(prefix))
os.system('mkdir -p %s/fordevelop'%(prefix))
os.system('mv %s.gene_name.notconsistent.list.txt  %s/fordevelop'%(prefix,prefix))
os.system('mv %s.gene_start.notconsistent.list.txt  %s/fordevelop'%(prefix,prefix))
os.system('mv %s.RefSeq_geneID.txt  %s/fordevelop'%(prefix,prefix))
print 'fastadic finished!!!'

os.system('cat /dev/null > %s.fasta'%(prefix))
g=open('%s.fasta'%(prefix),'r+')

# ////////////////////////////////////////////////////
# design started
os.system('cat /dev/null > %s.N20.fasta'%(prefix))
os.system('cat /dev/null > %s.N20NGG.fasta'%(prefix))
os.system('cat /dev/null > %s.gene_statistics.txt'%(prefix))
os.system('cat /dev/null > %s.sgRNA_statistics.txt'%(prefix))
m=open('%s.gene_statistics.txt'%(prefix),'r+')
n=open('%s.sgRNA_statistics.txt'%(prefix),'r+')
n.write('sgRNAID\tsgRNA_position_in_gene\tGCcontent\n')
m.write('gene_name\tgene_length\tsgRNA_number_in_gene\n')
m.close()
n.close()

for geneID in target_fasta_dic:
    print 'Design process for %s starts.  '%(geneID)
    gene_seq=target_fasta_dic[geneID]
    g.write('>%s\n%s\n'%(geneID,gene_seq))
    os.system('python find_N20_evaluate.py %s %s %s %s %s %s %s %s %s %s %s'%(geneID,gene_seq,genome,sgRNA_number,ORFcutoff,GCcontentMin,GCcontentMax,strand,off_threshold,spacer_length,PAM))
    os.system('cat %s_candidateN20.fasta >> %s.N20.fasta'%(geneID,prefix))
    os.system('cat %s_candidateN20NGG.fasta >> %s.N20NGG.fasta'%(geneID,prefix))
    os.system('cat %s_statistics.txt >> %s.gene_statistics.txt'%(geneID,prefix))
    os.system('cat %s_sgRNA_statistics.txt >> %s.sgRNA_statistics.txt'%(geneID,prefix))
    os.system('rm %s_candidateN20.fasta'%(geneID))
    os.system('rm %s_candidateN20NGG.fasta'%(geneID))
    os.system('rm %s_statistics.txt'%(geneID))
    os.system('rm %s_sgRNA_statistics.txt'%(geneID))
        
    
g.close()

geneLst=[]
positionLst=[]
with open('%s.gene_statistics.txt'%(prefix),'r') as f:
    for i,line in enumerate(f):
        if i!=0:
            line=line.split('\t')
            geneLst.append(int(line[2]))
geneLst.sort()
with open('%s.sgRNA_statistics.txt'%(prefix),'r') as f:
    for i,line in enumerate(f):
        if i!=0:
            line=line.split('\t')
            positionLst.append(float(line[1]))
positionLst.sort()

os.system('cat %s.N20.fasta > %s.N20.txt'%(prefix,prefix))
os.system('cat %s.N20NGG.fasta > %s.N20NGG.txt'%(prefix,prefix))
os.system('cat %s.fasta  > %s.txt'%(prefix,prefix))
os.system('mv %s.* %s'%(prefix,prefix))


plt.hist(geneLst,color='r')
plt.xlabel('sgRNA number per gene')
plt.ylabel('Gene number')
plt.title('The distrution of sgRNA number per gene')
plt.savefig('%s/The_distrution_of_sgRNA_number_per_gene.png'%(prefix),dpi=1000)
plt.clf()

plt.hist(positionLst,color='b')
plt.xlabel('sgRNA position(relative to CDS 5`)')
plt.ylabel('sgRNA number')
plt.title('The distrution of sgRNA position per gene')
plt.savefig('%s/The_distrution_of_sgRNA_position_per_gene.png'%(prefix),dpi=1000)

if negative=='yes':
    os.system('python negativeControl.py %s %s %s %s %s %s'%(genome,negative_number,spacer_length,GCcontentMin,GCcontentMax,prefix))
    os.system('mkdir -p %s/negative'%(prefix))
    os.system('cat %s_N20_NC_passed.fasta >%s/negative/%s_N20_NC_passed.txt'%(prefix,prefix,prefix))
    os.system('mv %s_N20_NC_passed.fasta  %s/negative'%(prefix,prefix))
