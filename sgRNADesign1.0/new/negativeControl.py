# this script is sued to construct the N20 library as negative control in the library
# principle:
# find N20 with defined structure:
# define a 2mer library with all the possible constrcut
# the content of 4mer should be given as a random distribution
# the GC content should be well controled. (0.25,0.75)
# Meanwhile, a search against genome is conducted by standalone seqmap tool to remove those with potential offtarget effect

# input: genome nucl seq in fasta format

# input: number of N20 negative control need to be designed

import os
import sys
import numpy as np

genome=sys.argv[1]
N20_number=int(sys.argv[2])
GCcontentMin=int(sys.argv[3])
GCcontentMax=int(sys.argv[4])
prefix=sys.argv[5]
output=sys.argv[6]
spacer_length=20
Twomer_list=['AA','AG','AT','AC','TT','TA','TC','TG','CC','CA','CG','CT','GA','GC','GG','GT']

Onemer_list=['A','T','C','G']
# obtain the GC content of N20
def GCcontent_extract(seq):
    seqLen=len(seq)
    GCnumber=0
    for nulc in seq:
        if nulc=='G' or nulc=='C':
            GCnumber+=1
    GContentPerc=100*float(GCnumber)/float(seqLen)
    return GContentPerc

# /////////////////////////////////////////////////////////////////////////////////////////////
# dic to store all the candidate N20NGG
# {'ID':seq, ...}
Candidate_N20={}

# put all the potential N20NGG into a fasta file and a dictionary with N20NGG as index
# candidateoff.fasta to store all the N20-N(abbr)G(A)G seq in fasta format
os.system('cat /dev/null > N20NGGofftarget.fasta')
sgRNA_ID_collection=[]

# //////////////////////////////
# get the initial N20 negative control list
(quotient,remainder)=divmod(spacer_length,2)
for two_mer in range(N20_number*4):
    random_N20=''
    for i in range(quotient):
        # get a random 2-mer for totally 10 times
        random_2mer=Twomer_list[np.random.randint(len(Twomer_list))]
        random_N20+=random_2mer
    if remainder==1:
        random_N20+=np.random.choice(Onemer_list)
    
    sgRNA_ID='NC_%d'%(two_mer)
    if sgRNA_ID not in sgRNA_ID_collection:
        sgRNA_ID_collection.append(sgRNA_ID)
    Candidate_N20[sgRNA_ID]=random_N20
    target1Nameline='>%s-Target1'%(sgRNA_ID)
    potentialOffTarget1=random_N20+'NGG'
    target2Nameline='>%s-Target2'%(sgRNA_ID)
    potentialOffTarget2=random_N20+'NAG'
    os.system('echo "%s\n%s\n%s\n%s\n" >> N20NGGofftarget.fasta'%(target1Nameline,potentialOffTarget1,target2Nameline,potentialOffTarget2))

# conduct seqmap to find the potential offtarget
os.system('./seqmap 5 N20NGGofftarget.fasta %s seqmapOut.temp.dat /output_all_matches'%(genome))

# ///////////////////////////////////////
# those gRNA without any hit would be selected

os.system('wc -l seqmapOut.temp.dat')

f=open('seqmapOut.temp.dat','r')
# store sgRNA with hits towards the genome detected by seqmap
sgRNA_withHit=[]

# to put the sgRNA passing the threshold into the doc
# to find the mismatch hit in genome for each sgRNA
for line in f:
    row=line.rstrip().split('\t')
    if row[0]!='trans_id':
        sgRNA_ID=row[3].split('-')[0]
        if sgRNA_ID not in sgRNA_withHit:
            sgRNA_withHit.append(sgRNA_ID)
f.close()

print 'number of sgRNA:    ',
print len(Candidate_N20)
print 'number of sgRNA_withHit:    ',
print len(sgRNA_withHit)

os.system('rm seqmapOut.temp.dat N20NGGofftarget.fasta')

#print sgRNA_withHit

# sgRNA that qulitified
number_record=0
g=open(output,'a')
tem=open("tem.txt",'w')
for sgRNA_ID in Candidate_N20:
#    print sgRNA_ID
    random_N20=Candidate_N20[sgRNA_ID]
    if sgRNA_ID not in sgRNA_withHit: # pass the quatification that no hit found in the seqmap search
        if GCcontent_extract(random_N20)<GCcontentMax and GCcontent_extract(random_N20)>GCcontentMin:  # GC content ok
             g.write('>%s\n%s\n'%(sgRNA_ID,random_N20))
             tem.write('>%s\n%s\n'%(sgRNA_ID,random_N20))
             number_record+=1
#             print number_record
    if number_record>=N20_number:
        break
g.close()
tem.close()
