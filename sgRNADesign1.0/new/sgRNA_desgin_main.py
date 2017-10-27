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
import ConfigParser

#The function processing the configure file
def configprocess(text):
    config=ConfigParser.ConfigParser()
    config.read(text)
    if len(config.sections())>1:
        print('Make a mistake in the configure file') and os._exit(0) 
    section=''.join(config.sections())
    varlist={option:config.get(section,option) for option in config.options(section)}
    return varlist
#input inspection 
def metricJudge(varlist):
    orfcutoff    =True if float(varlist['orfcutoff'])<1 else False
    sgrna_number =True if varlist['sgrna_number'].isdigit() else False
    gccontentmin =True if float(varlist['gccontentmin'])<100 and float(varlist['gccontentmin'])>0 else False
    gccontentmax =True if float(varlist['gccontentmax'])<100 and float(varlist['gccontentmax'])>0 else False
    off_threshold=True if float(varlist['off_threshold'])<1000 else False
    strand       =True if varlist['strand'] in ['nontemplate','template'] else False
    negative     =True if varlist['negative'] in ['Yes','No'] else False
    negative_number=True if varlist['negative_number'].isdigit() else False
    multiple     =True if varlist['multiple'] in ['Yes','No'] else False
    genomewide   =True if varlist['genomewide']  in ['Yes','No'] else False
    targetfasta  =True
    indexfile    =True
    genome       =True
    blastresult  =True
    prefix       =True
    for item in varlist:
        if not eval(item) :
            print('The varaint %s is wrong,please input the correct value!!!!'%(item))
            print(varlist[item])
            os._exit(1)
          
    
    
# processing the ptt or rnt file
def pttprocess(text):
    genePositionLst=[]
    genePositionDic={}
#    geneStrand={}
    warning='Sorry, the gene is neither in the template chain nor on the non-template chain'
    with open(text,'r') as pttfile:
        for record in pttfile.readlines()[3:]:
            items=record.rstrip().split('\t')
            synonym = items[5]  if '_' not in items[5] else items[5].replace('_','')  
            name =items[4]  if '_' not in items[4] else items[4].replace('_','')
            name =name  if ' ' not in name else name.replace(' ','') 
            name =name  if '-' not in name else name.replace('-','')
            start = np.where(items[1]=='+',items[0].split('..')[0],np.where(items[1]=='-',\
                             items[0].split('..')[1],warning))        
            genePositionLst.append(str(start))
            genePositionDic[str(start)]=name+synonym
    return genePositionLst,genePositionDic  

#judge the gene cluster that has a gene is similar to another cluster
def judge(TLst,OLst):
    flag=False
    for Lst in TLst:
        for item in OLst:
            if item in Lst:
                TLst[TLst.index(Lst)]=list(set(Lst+OLst))
                flag=True
                break
        if flag:
            break
    return flag
#find the similar gene and put them in a cluster 
def Lstprocess(recDic):
    record=[[item[0]]+item[1] for item in recDic.items()]
    group=[]# record the processed gene cluster
    combineLst=[] #output the list of gene cluster!
    for i in range(len(record)-1):
        if i not in group:
            combineLst.append(record[i])
            group.append(i)
            for j in range(i+1,len(record)):
                if j not in group:
                    flag=judge(combineLst,record[j])
                    if flag:
                        group.append(j)
        else:
            continue
    if len(record)-1 not in group:
        combineLst.append(record[len(record)-1])
    itemLst=[]
    for Lst in combineLst:
        for item in Lst:
            itemLst.append(item)
    return itemLst,combineLst

#process the blast file
def blastfileprocess(blastfile):
    recordLst=[]
    recordDic={}
    with open(blastfile)as f:
        for line in f.readlines():
            ID1=line.split('\t')[0].split(':')[1]
            ID2=line.split('\t')[1].split(':')[1]
            if ID1 not in recordLst:
                recordLst.append(ID1)
                recordDic[ID1]=[ID2]
            else:
                recordDic[ID1].append(ID2)
    MultiRecordDic={}
    for item in recordDic:
        if len(recordDic[item])>1:
            key=item.split('-')[0].split('c')[1] if 'c' in item else item.split('-')[0]
            value=[]
            for record in recordDic[item]:
                value.append(record.split('-')[0].split('c')[1] if 'c' in record else record.split('-')[0])
            MultiRecordDic[key]=value
    itemLst,combineLst=Lstprocess(MultiRecordDic)
    return itemLst,combineLst    

#Transform the position name to gene name
def nametransform(genePositionLst,genePositionDic,clusterLst,clusterGroupLst):
    clusterLst_trans=[]
    clusterGroupLst_trans=[]
    clusterLst_trans=[genePositionDic[position]  for position in clusterLst if position in genePositionLst ]
    for cluster in clusterGroupLst:
        sub=[genePositionDic[position] for position in cluster if position in genePositionLst ]
        clusterGroupLst_trans.append(sub)
    return clusterLst_trans,clusterGroupLst_trans

# construct a two-layer list with geneID where highly homologous gene clusters are grouped together as multi-member list, while others remain single-member list
# two layer list: if one gene has no highly homologous partner in the genome, append [gene], else, append [gene1, gene2, gene3, gene4, ...]        
def construct_two_layer(targetDic,clusterLst_trans,clusterGroupLst_trans):
    geneID_clusterLst=[]
    for geneID in targetDic:
        if geneID not in clusterLst_trans:
            geneID_clusterLst.append([geneID])
    return clusterGroupLst_trans+geneID_clusterLst

# start to check each entry of target fasta file
def tarFastaprocess(text,genePositionLst,genePositionDic):
    targetfastaDic={}
    with open(text,'r') as targetfastafile:
        for line in targetfastafile.readlines():
            if '>' in line:
                flag=False
                location_item=line[line.find(':')+1:line.find('-')] 
                start=location_item.replace('c','') if 'c' in location_item else location_item
                if start in genePositionLst:
                    geneID=genePositionDic[start] # to store the protein ID now
                    targetfastaDic[geneID]=''
                    flag=True
            elif flag:
                targetfastaDic[geneID]+=line.rstrip()
    
    return targetfastaDic    

#Design the sgRNA for Gene
def design(geneID,gene_seq,varlist):
    os.system('python find_N20_evaluate.py %s %s %s %s %s %s %s %s %s %s'%(geneID,gene_seq,\
              varlist['genome'],varlist['sgrna_number'],varlist['orfcutoff'],\
              varlist['gccontentmin'],varlist['gccontentmax'],varlist['strand'],\
              varlist['off_threshold'],varlist['multiple']))


#Plot the statistic of the sgRNA position information in one Gene and the sgRNA number information of the Gene
def resultplot(data,xlab,ylab,picturetitle,name,color):
    data.sort()    
    plt.hist(data,align='mid',color=color)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(picturetitle)
    plt.savefig(name,dpi=1000)
    plt.clf()

#Design the negative control sgRNA for the genome
def NegativeControl(varcontrol,output):
    os.system('python negativeControl.py %s %s %s %s %s %s'\
             %(varcontrol['genome'],varcontrol['negative_number'],varcontrol['gccontentmin'],\
               varcontrol['gccontentmax'],varcontrol['prefix'],output))
def nongenomic(varlist):
    targetfastaDic={}
    geneID_clusterLst=[]
    with open(varlist['targetfasta'],'r') as targetfastafile:
        for line in targetfastafile.readlines():
            if '>' in line:
                flag=False
                geneID=line.rstrip().split('>')[1]
                targetfastaDic[geneID]=''
                geneID_clusterLst.append([geneID])
                flag=True
            elif flag:
                targetfastaDic[geneID]+=line.rstrip()
    geneID_clusterLst=[[geneID] for geneID in targetfastaDic]
    return geneID_clusterLst,targetfastaDic

def singlecopy(varlist):
    genePositionLst,genePositionDic=pttprocess(varlist['indexfile'])
    target_fasta_dic=tarFastaprocess(varlist['targetfasta'],genePositionLst,genePositionDic)
    geneID_clusterLst=[[geneID] for geneID in target_fasta_dic]
    return geneID_clusterLst,target_fasta_dic
def muticopy(varlist):
    genePositionLst,genePositionDic=pttprocess(varlist['indexfile'])
    target_fasta_dic=tarFastaprocess(varlist['targetfasta'],genePositionLst,genePositionDic)
    os.system('python blast_filter.py %s %s %s'%(varlist['targetfasta'],varlist['targetfasta'],varlist['blastresult']))
    subset_Lst,subset_GroupLst=blastfileprocess('%s.trim'%(varlist['blastresult']))
    subset_transLst,subset_transGroupLst=nametransform(genePositionLst,genePositionDic,subset_Lst,subset_GroupLst)
    geneID_clusterLst=construct_two_layer(target_fasta_dic,subset_transLst,subset_transGroupLst)
    os.system('rm %s.trim'%(varlist['blastresult']))
    return geneID_clusterLst,target_fasta_dic
def processPositionfile(sgRNAPositionfile,sgRNA_position,writetype):
    w=open(sgRNA_position,writetype)
    with open(sgRNAPositionfile,'r') as f:
        for i,line in enumerate(f.readlines()[1:]):
            line=line.split('\t')
            text=line[0].split('_')[0]+'\t'+line[0]+'\t'+line[1]
            w.write('%s\n'%(text)) if i!=len(f.readlines())-1 else w.write('%s'%(text))
def processSgRNASequence(N20file,N20_library,writetype,NC):
    w=open(N20_library,writetype)
    if not NC:
        w.write('sgRNAID,sgRNAseq,gene')  
    with open(N20file,'r') as f:
        for i,line in enumerate(f.readlines()):
            line=line.strip()
            if line[0]=='>':
                text=line[1:]
            else:
                if NC:
                    text=text+','+line+','+'0'
                else:
                    text=text+','+line+','+text.split('_')[0]
                w.write('\n%s'%(text))
def main(configureFile):
    varlist= configprocess(configureFile)
    metricJudge(varlist)
    varlist['multiple']=True if varlist['multiple']=='Yes' else False
    varlist['genomewide']=True if varlist['genomewide']=='Yes' else False
    if varlist['multiple'] and varlist['genomewide']:
        geneID_clusterLst,target_fasta_dic=muticopy(varlist)
    elif not varlist['multiple'] and varlist['genomewide']:
        geneID_clusterLst,target_fasta_dic=singlecopy(varlist)
    elif varlist['multiple'] and not varlist['genomewide']:
        geneID_clusterLst,target_fasta_dic=nongenomic(varlist)
    else:
        geneID_clusterLst,target_fasta_dic=nongenomic(varlist)
    print 'fastadic finished!!!'
 
    with open('%s.gene_statistics.txt'%(varlist['prefix']),'w') as m:
        m.write('gene_name\tgene_length\tsgRNA_number_in_gene\n')
    with open('%s.sgRNA_statistics.txt'%(varlist['prefix']),'w') as n:
        n.write('sgRNAID\tsgRNA_position_in_gene\tGCcontent\n')
    os.system('cat /dev/null > %s.N20.fasta.txt'%(varlist['prefix']))
    os.system('cat /dev/null > %s.N20NGG.fasta.txt'%(varlist['prefix']))
    k=open('%s.cluster.txt'%(varlist['prefix']),'w')
    with open('%s.fasta.txt'%(varlist['prefix']),'w') as g:
        for geneCluster in geneID_clusterLst:
            geneID=geneCluster[0]
            group=','.join(geneCluster)
            gene_seq=''
            for gene in geneCluster:
                gene_seq+=target_fasta_dic[gene]+'-'
                g.write('>%s\n%s\n'%(gene,target_fasta_dic[gene]))
#                os._exit(1)
            gene_seq=gene_seq[:-1]
#            os._exit(1)

            k.write('%s\t%s\n'%(geneID,group)) 
            print 'Design process for %s starts.  '%(geneID)
            print 'Copy number: %d'%(len(geneCluster))
            print 'Gene length: %s'%(float(len(gene_seq))/float(len(geneCluster)))
            design(geneID,gene_seq,varlist)
            os.system('cat %s_candidateN20.fasta >> %s.N20.fasta.txt'%(geneID,varlist['prefix']))
            os.system('cat %s_candidateN20NGG.fasta >> %s.N20NGG.fasta.txt'%(geneID,varlist['prefix']))
            os.system('cat %s_statistics.txt >> %s.gene_statistics.txt'%(geneID,varlist['prefix']))
            os.system('cat %s_sgRNA_statistics.txt >> %s.sgRNA_statistics.txt'%(geneID,varlist['prefix']))
            os.system('rm %s_*'%(geneID))
    k.close()
    sgRNAfile=os.popen('cut -f 3 %s.gene_statistics.txt'%(varlist['prefix'])).readlines()[1:]
    sgRNAfile=[int(sgRNAnumer.rstrip()) for sgRNAnumer in sgRNAfile]
    positionfile=os.popen('cut -f 2 %s.sgRNA_statistics.txt'%(varlist['prefix'])).readlines()[1:]
    positionfile=[float(sgRNAposition.rstrip()) for sgRNAposition in positionfile]
    resultplot(sgRNAfile,'sgRNA number per gene','Gene number','The distrution of sgRNA number per gene','The_distrution_of_sgRNA_number_per_gene.png','r')
    resultplot(positionfile,'sgRNA position(relative to CDS 5`)','sgRNA number','The distrution of sgRNA position per gene','The_distrution_of_sgRNA_position_per_gene.png','b')
    os.system('mkdir %s'%(varlist['prefix']))
    os.system('mv The_distrution_of* %s'%(varlist['prefix']))
    processSgRNASequence('%s.N20.fasta.txt'%(varlist['prefix']),'N20_library.csv','w',False)
    processPositionfile('%s.sgRNA_statistics.txt'%(varlist['prefix']),'sgRNA_position.txt','w')
    if varlist['negative']=='Yes':
        NegativeControl(varlist,'%s.N20.fasta.txt'%(varlist['prefix']))
        processSgRNASequence('tem.txt','N20_library.csv','a',True)
        w=open('sgRNA_position.txt','a')
        with open('tem.txt','r') as f:
            for i,line in enumerate(f.readlines()):
                line=line.strip()
                if '>' in line:
                    text='0'+'\t'+line[1:]+'\t'+'0'
                    w.write('%s\n'%(text)) if i!=len(f.readlines())-2 else w.write('%s'%(text))
        os.system("rm tem.txt")
    os.system('mv N20_library.csv %s'%(varlist['prefix']))
    os.system('mv sgRNA_position.txt %s'%(varlist['prefix']))
    os.system('mv %s.* %s'%(varlist['prefix'],varlist['prefix']))


if __name__ =='__main__':
    main(sys.argv[1])
    
