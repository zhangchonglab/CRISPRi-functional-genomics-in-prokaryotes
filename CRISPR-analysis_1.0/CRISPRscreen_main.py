# this script is used as the main function for the CRISPR screening pipeline
# based on this script, we can choose to do the pretreatment, sgRNA level statistics calculation, gene level statistics calculation and operon level statistics calculation, respectively.

import os
import sys
config_file=sys.argv[1]
config={}
with open (config_file,'r') as f:
    for line in f:
        line=line.strip()
        if len(line)==0:
            continue
        if line[0]=='#':
            continue
        row=line.split('\t')
        config[row[0]]=[]
        if len(row)==1:
            config[row[0]]=''
        elif len(row)==2:
            config[row[0]]=row[1]
        else:
            config[row[0]]=row[1:]
#//////////////////////////////////////////////////////////////
#process the configure file 
# prefix for naming of all files
if 'prefix' not in config or config['prefix']=='':
    prefix='Test'
    print('the default value of prefix is Test!')
else:
    prefix=config['prefix']
#read the path of fastq file 
if 'fastqpath' not in config:
    print('please input the path of fastq file' )
    sys.exit()
else:
    if config['fastqpath'][-1]!='/':
        print('please input the correct path')
        sys.exit()
    else:
        fastqpath=config['fastqpath']
#read the fastq file and corresponding label
for i in range(len(config['fastq'])):
    
    config['fastq'][i]=os.path.expanduser(fastqpath)+config['fastq'][i]
#fastq=','.join(config['fastq']).replace(',',' ')
fastq=','.join(config['fastq'])
#read the corresponding label of fastq file
if 'sample-label' not in config or config['sample-label']=='':
    print('please input the sample label')
    sys.exit()
else:
    sample_label=','.join(config['sample-label'])

#read the primer sequence
if 'forward_prefixseq' not in config or config['forward_prefixseq']=='':
    print('Please input the forward_prefixseq')
    sys.exit()
if 'forward_suffixseq' not in config or config['forward_suffixseq']=='':
    print('Please input the forward_suffixseq')
    sys.exit()
forward_prefixseq=config['forward_prefixseq']
forward_suffixseq=config['forward_suffixseq']

#specify the length of the sgRNAs (without PAM sequence)
if 'sgrna-len' not in config or config['sgrna-len']=='':
    sgrna_len=20
    print('the default value of sgrna_len is 20!')
else:
    sgrna_len=config['sgrna-len']

#Path to library design file (csv format, columns: id, sequence, gene)
if 'list-seq' not in config or config['list-seq']=='':
    print('please input the sgRNA library file!')
    sys.exit()
else:
    list_seq=config['list-seq']
#experiment_configure=sys.argv[1]
if 'experiment_configure' not in config or config['experiment_configure']=='':
    print('please input the experiment design file!')
    sys.exit()
else:
    experiment_configure=config['experiment_configure']


# method to set the control sgRNA set: 'all or NC'
if 'control_setting' not in config or config['control_setting']=='' :
    control_setting='all'
    print('the default value of control_setting is all!')
else:
    if config['control_setting']=='all' or config['control_setting']=='NC':
        control_setting=config['control_setting']
    else:
        print('Please input the correct method selection!')
        sys.exit()

# FDR threshold for hit gene calling
if 'FDR_threshold' not in config or config['FDR_threshold']=='':
    FDR_threshold=0.05
    print('The default value of FDR_threshold is 0.05!')
else:
    FDR_threshold=config['FDR_threshold']
    if float(FDR_threshold)>1:
        print('The value of FDR_threshold must be less than 1!')
        sys.exit()


#number of reads in initial library that below this threshold won't be incoperated in to the following analysis pipleine 
if 'ReadsThreshold' not in config or config['ReadsThreshold']=='':
    ReadsThreshold=100
    print('the default value of ReadsThreshold is 100!')
else:
    ReadsThreshold=config['ReadsThreshold']

# method for hit gene calling
if 'hit_gene_calling' not in config or config['hit_gene_calling']=='':
    hit_gene_calling='all'
    print('The default value of hit_gene_calling is all!')
else:
    if config['hit_gene_calling']=='all' or config['hit_gene_calling']=='position':
        hit_gene_calling=config['hit_gene_calling']
    else:
        print('Please input the correct method for hit gene calling!')
        sys.exit()
#flat file of gene sgRNA position information
# sgRNA gene ORF_position(%)
# tab deliminated
gene_sgRNA_position=config['gene_sgRNA_position']

#Operon_gene_List
if 'Operon_gene_List'not in config or config['Operon_gene_List']=='':
    Operon_gene_List=''
else:
    Operon_gene_List=config['Operon_gene_List']


#quasi gene number for one condition
if 'QuasiGene_number' not in config or config['QuasiGene_number']=='':
    QuasiGene_number=1000
#    print('the default value of QuasiGene_number is 1000!')
else:
    QuasiGene_number=config['QuasiGene_number']

# minimal sgRNA number for p Value calculation in subsampling
if 'min_sgRNAnumber' not in config or config['min_sgRNAnumber']=='':
    min_sgRNAnumber=5
#    print('the default value of min_sgRNAnumber is 5!')
else:
    min_sgRNAnumber=config['min_sgRNAnumber']
#minimal number and maxium number os sgRNAs for one quasi gene
if 'sgRNAperQuasiGene' not in config or config['sgRNAperQuasiGene']=='':
    sgRNAperQuasiGene='1,40'
#    print('the default minimal and maxium value of sgRNAperQuasiGene is 1 and 40!')
else:
    sgRNAperQuasiGene=config['sgRNAperQuasiGene'][0]+','+config['sgRNAperQuasiGene'][1]




# ////////////////////////////////////////////////////////////////////
# normalize the maped sgRNA number
outputname='%s_results/%s_count/%s'%(prefix,prefix,prefix)
os.system('python CRISPRscreen_normalize.py --list_seq %s --sample_label %s --output_prefix %s --prefix_nucl %s --suffix_nucl %s --variable_region_len %s --unmapped-to-file --fastq %s  2>error.log'%(list_seq,sample_label,outputname,forward_prefixseq,forward_suffixseq,sgrna_len,fastq))
print('normalize finished')
# ////////////////////////////////////////////////////////////////////
# sgRNA level statistics calculation
normalizedData='%s_results/%s_count/%s.normalizeCount.txt'%(prefix,prefix,prefix)
os.system('python CRISPRscreen_sgRNA.py %s %s %s %s %s 2>>error.log'%(experiment_configure,normalizedData,control_setting,prefix,ReadsThreshold))
geneDic='%s_results/%s_Processed_geneDic.pickle'%(prefix,prefix)
geneZDic='%s_results/%s_Processed_geneZDic.pickle'%(prefix,prefix)
geneLst='%s_results/%s_processed_geneLst.pickle'%(prefix,prefix)
sgRNALst='%s_results/%s_processed_sgRNALst.pickle'%(prefix,prefix)
stressed_conditionLst='%s_results/%s_stressed_conditionLst.pickle'%(prefix,prefix)
print ('sgRNA statistics calculation finalized')
# ////////////////////////////////////////////////////////////////////
# gene level statistics calculation
os.system('python CRISPRscreen_gene.py %s %s %s %s %s %s %s %s %s %s %s %s %s 2>>error.log'%(geneDic,geneZDic,geneLst,sgRNALst,stressed_conditionLst,gene_sgRNA_position,control_setting,hit_gene_calling,FDR_threshold,prefix,sgRNAperQuasiGene,QuasiGene_number,min_sgRNAnumber))
Finish_result_geneDic='%s_results/%s_Result_geneDic.pickle'%(prefix,prefix)  #
print ('gene statistics calculation finalized')
#sys.exit()

# ////////////////////////////////////////////////////////////////////
# operon level statistics calculation
if len(Operon_gene_List)!=0:
    os.system('python CRISPRscreen_Operon.py %s %s %s 2>>error.log'%(Finish_result_geneDic,Operon_gene_List,prefix))

os.system('rm %s_results/*.pickle'%(prefix))
