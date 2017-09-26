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

# processing the ptt or rnt file
def pttprocess(text):
    genePositionLst=[]
    genePositionDic={}
    warning='Sorry, the gene is neither in the template chain nor on the non-template chain'
    with open(text,'r') as pttfile:
        for record in pttfile.readlines()[3:]:
            items=record.rstrip().split('\t')
            synonym = items[5]  if '_' not in items[5] else items[5].replace('_','')   
            start = np.where(items[1]=='+',items[0].split('..')[0],np.where(items[1]=='-',\
                             items[0].split('..')[1],warning))        
            genePositionLst.append(str(start))
            genePositionDic[str(start)]=synonym 
    return genePositionLst,genePositionDic  

#process the blast file
def blastfileprocess(blastfile):
    iden90_qC90_hC90_GroupLst=[]
    iden90_qC90_hC90_Lst=[]
    with open(blastfile,'r') as g:
        for line in g:
            row=line.rstrip().split('\t')
            gene1=row[0][row[0].find(':'):row[0].find('-')] if 'c' not in row[0] \
             else row[0][row[0].find(':'):row[0].find('-')].replace('c','')
            gene2=row[1][row[1].find(':'):row[1].find('-')] if 'c' not in row[1]\
             else row[1][row[1].find(':'):row[1].find('-')].replace('c','')
            if gene1!=gene2:
                if gene1 not in iden90_qC90_hC90_Lst and gene2 not in iden90_qC90_hC90_Lst:
                    iden90_qC90_hC90_Lst.append(gene1)
                    iden90_qC90_hC90_Lst.append(gene2)
                    iden90_qC90_hC90_GroupLst.append([gene1,gene2])
                elif gene1 in iden90_qC90_hC90_Lst and gene2 not in iden90_qC90_hC90_Lst:
                    iden90_qC90_hC90_Lst.append(gene2)
                    for geneCluster in iden90_qC90_hC90_GroupLst:
                        if gene1 in geneCluster:
                            geneCluster.append(gene2)
                elif gene2 in iden90_qC90_hC90_Lst and gene1 not in iden90_qC90_hC90_Lst:
                    iden90_qC90_hC90_Lst.append(gene1)
                    for geneCluster in iden90_qC90_hC90_GroupLst:
                        if gene2 in geneCluster:
                            geneCluster.append(gene1)
    os.remove(blastfile)
    return iden90_qC90_hC90_Lst,iden90_qC90_hC90_GroupLst

#Transform the position name to gene name
def nametransform(genePositionLst,genePositionDic,iden90_qC90_hC90_Lst,iden90_qC90_hC90_GroupLst):
    iden90_qC90_hC90_Lst_trans=[]
    iden90_qC90_hC90_GroupLst_trans=[]
    iden90_qC90_hC90_Lst_trans=[genePositionDic[position]  for position in iden90_qC90_hC90_Lst if position in genePositionLst ]
    for cluster in iden90_qC90_hC90_GroupLst:
        sub=[genePositionDic[position] for position in cluster if position in genePositionLst ]
        den90_qC90_hC90_GroupLst_trans.append(sub)
    return iden90_qC90_hC90_Lst_trans,iden90_qC90_hC90_GroupLst_trans

# construct a two-layer list with geneID where highly homologous gene clusters are grouped together as multi-member list, while others remain single-member list
# two layer list: if one gene has no highly homologous partner in the genome, append [gene], else, append [gene1, gene2, gene3, gene4, ...]        
def construct_two_layer(targetDic,iden90_qC90_hC90_Lst_trans,iden90_qC90_hC90_GroupLst_trans):
    geneID_clusterLst=[]
    for geneID in targetDic:
        if geneID not in iden90_qC90_hC90_Lst_trans:
            geneID_clusterLst.append([geneID])
    return geneID_clusterLst+iden90_qC90_hC90_GroupLst_trans


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
    os.system('python find_N20_evaluate.py %s %s %s %s %s %s %s %s %s'%(geneID,gene_seq,\
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
def NegativeControl(varcontrol):
    os.system('python negativeControl.py %s %s %s %s %s '\
             %(varcontrol['genome'],varcontrol['negative_number'],varcontrol['gccontentmin'],\
               varcontrol['gccontentmax'],varcontrol['prefix']))
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
    geneID_clusterLst=[[geneID] for geneID in target_fasta_dic]
    return geneID_clusterLst,target_fasta_dic

def singlecopy(varlist):
    genePositionLst,genePositionDic=pttprocess(varlist['indexfile'])
    target_fasta_dic=tarFastaprocess(varlist['targetfasta'],genePositionLst,genePositionDic)
    geneID_clusterLst=[[geneID] for geneID in target_fasta_dic]
    return geneID_clusterLst,target_fasta_dic
def muticopy(varlist):
    genePositionLst,genePositionDic=pttprocess(varlist['indexfile'])
    target_fasta_dic=tarFastaprocess(varlist['targetfasta'],genePositionLst,genePositionDic)
    os.system('python blast_filter.py %s %s %s'%(varlist['targetfasta'],varlist['targetfasta'],varlist['blastresult']))
    os.system('cut -f 1 %s.trim.iden90.queCo90.hitCo90 |sort |uniq -d > temp.txt'%(varlist['blastresult']))
    os.system('fgrep -wf temp.txt %s.trim.iden90.queCo90.hitCo90 > %s.trim.iden90.queCo90.hitCo90.duplicate'%(varlist['blastresult'],varlist['blastresult']))
    os.remove('temp.txt')
    os.remove('%s.trim.iden90.queCo90.hitCo90'%(varlist['blastresult']))
    subset_Lst,subset_GroupLst=blastfileprocess('blastresult.trim.iden90.queCo90.hitCo90.duplicate')
    subset_transLst,subset_transGroupLst=nametransform(genePositionLst,genePositionDic,subset_Lst,subset_GroupLst)
    geneID_clusterLst=construct_two_layer(target_fasta_dic,subset_transLst,subset_transGroupLst)
    return geneID_clusterLst,target_fasta_dic

def main(configureFile):
    varlist= configprocess(configureFile)
    if varlist['multiple'] and varlist['nongenome']:
        print('The method is single copy')
        geneID_clusterLst,target_fasta_dic=nongenomic(varlist)
    elif varlist['multiple'] and not varlist['nongenome']:
        geneID_clusterLst,target_fasta_dic=muticopy(varlist) 
    elif not varlist['multiple'] and varlist['nongenome']:
        geneID_clusterLst,target_fasta_dic=nongenomic(varlist)
    else:
        geneID_clusterLst,target_fasta_dic=singlecopy(varlist)
    print 'fastadic finished!!!'
    with open('%s.gene_statistics.txt'%(varlist['prefix']),'w') as m:
        m.write('gene_name\tgene_length\tsgRNA_number_in_gene\n')
    with open('%s.sgRNA_statistics.txt'%(varlist['prefix']),'w') as n:
        n.write('sgRNAID\tsgRNA_position_in_gene\tGCcontent\n')
    os.system('cat /dev/null > %s.N20.fasta'%(varlist['prefix']))
    os.system('cat /dev/null > %s.N20NGG.fasta'%(varlist['prefix']))
    k=open('%s.cluster.txt'%(varlist['prefix']),'w')
    with open('%s.fasta'%(varlist['prefix']),'w') as g:
        for geneCluster in geneID_clusterLst:
            geneID=geneCluster[0]
            group=','.join(geneCluster)
            gene_seq=''
            for gene in geneCluster:
                gene_seq+=target_fasta_dic[gene]+'-'
                g.write('>%s\n%s\n'%(gene,target_fasta_dic[gene]))
            gene_seq=gene_seq[:-1]
            k.write('%s\t%s\n'%(geneID,group)) 
            print 'Design process for %s starts.  '%(geneID)
            print 'Copy number: %d'%(len(geneCluster))
            print 'Gene length: %s'%(float(len(gene_seq))/float(len(geneCluster)))
            design(geneID,gene_seq,varlist)
            os.system('cat %s_candidateN20.fasta >> %s.N20.fasta'%(geneID,varlist['prefix']))
            os.system('cat %s_candidateN20NGG.fasta >> %s.N20NGG.fasta'%(geneID,varlist['prefix']))
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
    os.system('cat %s.N20.fasta > %s.N20.txt'%(varlist['prefix'],varlist['prefix']))
    os.system('cat %s.N20NGG.fasta > %s.N20NGG.txt'%(varlist['prefix'],varlist['prefix']))
    os.system('cat %s.fasta  > %s.txt'%(varlist['prefix'],varlist['prefix']))
    os.system('mkdir %s'%(varlist['prefix']))
    os.system('mv %s.* %s'%(varlist['prefix'],varlist['prefix']))
    os.system('mv The_distrution_of* %s'%(varlist['prefix']))
    if varlist['negative']=='yes':
        NegativeControl(varlist)
        os.system('mkdir -p %s/negative'%(varlist['prefix']))
        os.system('cat %s_N20_NC_passed.fasta >%s/negative/%s_N20_NC_passed.txt'%(varlist['prefix'],varlist['prefix'],varlist['prefix']))
        os.system('mv %s_N20_NC_passed.fasta  %s/negative'%(varlist['prefix'],varlist['prefix']))



if __name__ =='__main__':
    main(sys.argv[1])
    
