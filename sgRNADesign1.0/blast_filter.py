# this script is used to extract the desired hits from a blast result based on several thresholds:
# identity, coverage for query, coverage for hit;

# first of all get the query and database fastafile, respectively
import os
import sys

query_file=sys.argv[1]
database_file=sys.argv[2]
blast_result_file=sys.argv[3]
# first of all, ask the user to define the threshold that they want to use:
# identity
identity_thre=True
identity=0.95
print 'Identity: %s'%str(identity)
# coverage for query
queryCoverage_thre=True
queryCoverage=0.95
print 'queryCoverage: %s'%str(queryCoverage)
# coverage for hit
hitCoverage_thre=True
hitCoverage=0.95
print 'hitCoverage: %s'%str(hitCoverage)

# now construct the query and hit length dic with data structure like this
# {queryID:float,queryID:float}
# the same as hit
# query_length_dic and hit_length_dic, respectively
query_length_dic={}
proteinID='' # to store the protein ID now
f=open(query_file,'r')
for line in f:
    if line[0]=='>':
        proteinID=line[1:(-1)].split(' ')[0]
        query_length_dic[proteinID]=0.0
    else:
        query_length_dic[proteinID]+=float(len(line.rstrip()))
f.close()        

hit_length_dic={}
proteinID='' # to store the protein ID now
f=open(database_file,'r')
for line in f:
    if line[0]=='>':
        proteinID=line[1:(-1)].split(' ')[0]
        hit_length_dic[proteinID]=0.0
    else:
        hit_length_dic[proteinID]+=float(len(line.rstrip()))
f.close()        

# now it is time for filtering
f=open(blast_result_file,'r')
g=open('%s.trim'%(blast_result_file),'w')
keys=['queryId', 'subjectId', 'percIdentity', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
#keys=['queryId', 'subjectId', 'percIdentity', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
blast_result=[]
for line in f:
    row=dict.fromkeys(keys)
    elements=line.rstrip().split('\t')
    for i, elem in enumerate(elements):
        row[keys[i]]=elem
    flag=True  # judge whether this result pass or not
    if identity_thre:
        if float(row['percIdentity'])/100.0<identity:
            flag=False
    if queryCoverage_thre:
        result_queryCoverage=float(row['alnLength'])/float(query_length_dic[row['queryId']])
        if result_queryCoverage<queryCoverage:
            flag=False
    if hitCoverage_thre:
        result_hitCoverage=float(row['alnLength'])/float(hit_length_dic[row['subjectId']])
        if result_hitCoverage<hitCoverage:
            flag=False
    if flag:
        g.write(line)
f.close()
g.close()







# trim the result and 
