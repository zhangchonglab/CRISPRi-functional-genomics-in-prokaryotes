from __future__ import print_function
import sys
import re
import os
import argparse
import math
import logging
import string
import subprocess
import matplotlib.pyplot as plt
import numpy as np


def CRISPRsgRNA_parseargs():
    """
    Parse arguments. Only used when CRISPRsgRNA.py is executed directly.
    """
    parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
    parser.add_argument('-l','--list_seq',required=True,help='A file containing the list of sgRNA names, their sequences and the gene in which sgRNA is . Support file format: csv and txt.')
    parser.add_argument('--sample_label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided. Default "sample1,sample2,...".')
    parser.add_argument('-n','--output_prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    parser.add_argument('--prefix_nucl',default='A',help='nucleotide sequence upstream the varaible region. Default A')
    parser.add_argument('--suffix_nucl',default='A',help='nucleotide sequence downstream the varaible region. Default A')
    parser.add_argument('--variable_region_len',type=int,default=20,help='Length of the sgRNA. Default 20')
    parser.add_argument('--unmapped-to-file',action='store_true',help='Save unmapped reads to file.')
    parser.add_argument('--fastq',default='',help='Sample fastq files, separated by comma (,). For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq,sample2_replicate1.fastq,sample2_replicate2.fastq" ')
    args=parser.parse_args()
    return args


def CRISPRsgRNA_revcomp(x):
    '''
    Reverse complement
    '''
    return x.translate(string.maketrans("ACGT","TGCA"))[::-1]

class VisualRCount:
    '''
    Class for generating reports of count command
    '''
    outprefix='sample1'

    # internal variable, for R file
    outrfh=None;  # file handle for R file

    # for Rnw file
    rnwtemplatestr=''
    outrnwfh=None
    outrnwstring=''
    # for statistics of coutns
    fastqfile=[]; # fastq files
    fastqlabels=[]; # fastq labels
    reads=[]; # read counts
    mappedreads=[]; # # mapped reads
    unmapduetosynthesis=[]; # # unmapped reads due to the synthesis errors
    unmapunknown=[]; # # unmapped reads with unknown sources
    totalsgrnas=[]; # # sgRNAs (in the library)
    zerocounts=[]; # # 0-count sgRNAs
    gini=[];  # gini index
    '''
    Member functions
    '''
    def setPrefix(self,prefix):
        '''
        Set up proper prefix
        '''
        (file_dir,file_base)=os.path.split(prefix)
        file_base_dot=re.sub('\.','_',file_base)
        self.outprefix=os.path.join(file_dir,file_base_dot)

    def writeCountSummary(self):
        '''
        Write statistics from gene summary file to buffer
        '''
        # insert string
        insertstr=''
        import textwrap
        fastqwrap=[' '.join(textwrap.wrap(x,50)) for x in self.fastqfile]
        insertstr+='filelist=c(' + ','.join(['"'+x+'"' for x in fastqwrap])  +');\n'
        insertstr+='labellist=c('+ ','.join(['"'+x+'"' for x in self.fastqlabels]) +');\n'
        insertstr+='reads=c('+','.join([str(x) for x in self.reads])+');\n'
        insertstr+='mappedreads=c('+','.join([str(x) for x in self.mappedreads])+');\n'
        insertstr+='unmap due to syn error=c('+','.join([str(x) for x in self.unmapduetosynthesis])+');\n'
        insertstr+='unmap due to unknown source=c('+','.join([str(x) for x in self.unmapunknown])+');\n'
        insertstr+='totalsgrnas=c('+','.join([str(x) for x in self.totalsgrnas])+');\n'
        insertstr+='zerocounts=c('+','.join([str(x) for x in self.zerocounts])+');\n'
        insertstr+='giniindex=c('+','.join([str(x) for x in self.gini])+');\n'
        #
        nwktowrite=re.sub('#__COUNT_SUMMARY_STAT__',insertstr,self.outrnwstring)
        # file names as list, instead of tables; disabled currently
        insertstr=''
        insertstr+='The fastq files processed are listed as follows.\n\n'
        insertstr+=r"\\"+"begin{enumerate}\n"
        for fq in self.fastqfile:
            insertstr+=r"\\"+"item "+fq+"\n"
        insertstr+=r"\\"+"end{enumerate}\n"
        # nwktowrite=re.sub('%__FILE_SUMMARY__',insertstr,nwktowrite)
        self.outrnwstring=nwktowrite

    def writeCountSummaryToTxt(self,txtfile):
        '''
        A stand-alone function to write the count summary to txt file
        fastqfile, fastqlabels, reads, mappedreads, and zerocounts must be set up
        '''
        ofstr=open(txtfile,'w')
        nsp=len(self.fastqfile)
        print('\t'.join(['File','Label','Reads','Mapped','Synerror','Unknown','Percentage','TotalsgRNAs',  'Zerocounts','GiniIndex']),file=ofstr)
        for i in range(nsp):
            print('\t'.join([self.fastqfile[i],self.fastqlabels[i],str(self.reads[i]),str(self.mappedreads[i]),str(self.unmapduetosynthesis[i]),str(self.unmapunknown[i]),"{:.4g}".format(self.mappedreads[i]*1.0/self.reads[i]),str(self.zerocounts[i]),"{:.4g}".format(self.gini[i])]),file=ofstr)
        #name.append('\t'.join([self.fastqfile[i])
        #gini.append("{:.4g}".format(self.gini[i])]))
        #per.append("{:.4g}".format(self.mappedreads[i]*1.0/self.reads[i]))
        ofstr.close()


def CRISPRsgRNA_checkargs(args):
    """
    Check args
    """
    if args.sample_label!='':
        nlabel=args.sample_label.split(',')
        nfq=args.fastq.split(',')
        if len(nlabel)!=len(nfq):
            logging.error('The number of labels ('+str(nlabel)+') must be equal to the number of fastq files provided.')
            sys.exit(-1)
    return 0

def CRISPRsgRNA_checklists(args):
    """
    Read sgRNA library file in csv or txt file
    format: sgRNAid  seq geneid
    Return
    sgRNA_geneDic,sgRNADic, sgRNALst
    ({sgRNAid,geneid}, {seq:id}, [id])
    """
    sgRNADic={}
    ctab={}
    sgRNA_geneDic={}
    n=0
    sgRNALst=[]
    for i,line in enumerate(open(args.list_seq)):
        
        field=line.strip().split(',')
        if i==0:
            continue
        n+=1
        if field[0] in sgRNADic:
            logging.warning('Duplicated sgRNA label '+field[0]+' in line '+str(n)+'. Skip this record.')
            continue
        if len(field)<3:
            logging.warning('Not enough field in line '+str(n)+'. Skip this record.')
            continue
        sgRNAid=field[0]
        sgRNAseq=field[1].upper()
        geneid=field[2]
        if sgRNAseq in sgRNALst:
            logging.warning('Duplicated sgRNA sequence '+field[1]+' in line '+str(n)+'. Skip this record.')
            continue
        sgRNADic[sgRNAseq]=sgRNAid
        sgRNALst.append(sgRNAid)
        sgRNA_geneDic[sgRNAid]=[geneid]
    logging.info('Loading '+str(len(sgRNALst))+' predefined sgRNAs.')
    return (sgRNA_geneDic,sgRNADic, sgRNALst)


#stastic partion function !!!!

def check_read(read,prefix,suffix,variableL,libDic):
    '''
    check whether the read contains member of the sgRNA library
    Parameters
    ----------
    filename
    read string
    prefix
    prefix seq upstream variable region
    suffix
    suffix seq downstream variable region
    variableL
    length of variable region
    libDic
    {sequence:id} dictionary
    Return value
    -----------
    id of identified library seq of the read
    '' if not detected
    '''
    pre_table=[n for n in xrange(len(read)) if read.find(prefix, n) == n]
    suf_table=[n for n in xrange(len(read)) if read.find(suffix, n) == n]
    sgRNAid='unknown'
    for p in [(x,y) for x in pre_table for y in suf_table]:
        if (p[1]-p[0]-len(prefix))>0.95*variableL and (p[1]-p[0]-len(prefix))<1.05*variableL:
            sgRNAid='synthesis error'
        if (p[1]-p[0]-len(prefix))==variableL:
            if read[(p[1]-variableL):p[1]] in libDic:
                sgRNAid=libDic[read[(p[1]-variableL):p[1]]]
                return sgRNAid
    rc_read=CRISPRsgRNA_revcomp(read)
    rc_pre_table=[n for n in xrange(len(rc_read)) if rc_read.find(prefix, n) == n]
    rc_suf_table=[n for n in xrange(len(rc_read)) if rc_read.find(suffix, n) == n]
    for p in [(x,y) for x in rc_pre_table for y in rc_suf_table]:
        if (p[1]-p[0]-len(prefix))>0.95*variableL and (p[1]-p[0]-len(prefix))<1.05*variableL:
            sgRNAid='synthesis error'
        if (p[1]-p[0]-len(prefix))==variableL:
            if rc_read[(p[1]-variableL):p[1]] in libDic:
                sgRNAid=libDic[rc_read[(p[1]-variableL):p[1]]]
                return sgRNAid
    return sgRNAid

#stastic partion function !!!!
def gettotalnormfactor(ctable):
    """
        Get the factor by total normalization
    """
    n=len(ctable[list(ctable.keys())[0]]) # samples
    m=len(ctable) # sgRNAs
    # calculate the sum
    sumsample=[0]*n
    for (k,v) in ctable.items():
        sumsample=[sumsample[i]+v[i] for i in range(n)]
    # normalizing factor
    avgsample=sum(sumsample)/float(n)
    samplefactor=[avgsample/k for k in sumsample]
    return samplefactor

def normalizeCount(ctable):
    """
        Central function for normalizing read counts
        Parameters:
        --------------
        ctable
        A dictionary of read counts: {sgrnaID:[count0,count1,...]}
        method
        Normalization methods: total
        Return value:
        --------------
        {sgRNA:[read counts]} if returnfactor == False, or [size_factor] if returnfactor == True
        
        By default, for higher read depths, the size factor is <1. If reversefactor is set, the factor is >1 (or 1/factor)
        """
    if len(ctable)==0:
        return ctable.copy()
    n=len(ctable[list(ctable.keys())[0]]) # samples
    m=len(ctable) # sgRNAs
    # calculate the total factor
    samplefactor=gettotalnormfactor(ctable)
    logging.debug('Initial (total) size factor: '+' '.join([str(x) for x in samplefactor]))
    logging.info('Final size factor: '+' '.join([str(x) for x in samplefactor]))
    ntable={ k: [ samplefactor[i]*v[i] for i in range(n)] for (k,v) in ctable.items()}
    return ntable


def CRISPRsgRNA_gini(x):
    '''
    Return the Gini index of an array
    Calculation is based on http://en.wikipedia.org/wiki/Gini_coefficient
    '''
    xs=sorted(x)
    n=len(xs)
    gssum=sum([ (i+1.0)*xs[i] for i in range(n)])
    ysum=sum(xs)
    if ysum==0.0:
        ysum=1.0
    gs=1.0-2.0*(n-gssum/ysum)/(n-1)
    return gs


def CRISPRsgRNA_mergedict(dict0,dict1):
    '''
    Merge all items in dict1 to dict0.
    dict0: {id:[count,count, ...], ...} the sequence of the count in list follows library in the args
    dict1: {id:count, id:count, ...}
    '''
    nsample=0
    if len(dict0)>0:
        nsample=len(dict0[list(dict0.keys())[0]])
    for (k,v) in dict0.items():
        if k in dict1:
            v+=[dict1[k]]
        else:
            v+=[0]
    for (k,v) in dict1.items():
        if k not in dict0:
            if nsample>0:
                dict0[k]=[0]*nsample
            else:
                dict0[k]=[]
            dict0[k]+=[v]
    # return dict0


def CRISPRsgRNA_processonefile(filename,args,mapDic,unmapDic,sgRNADic,datastat,prefix,suffix):
    '''
    Go through one fastq file
    Parameters
    ----------
    filename
    Fastq filename to be sequence
    args
    Arguments
    mapDic:
    A dictionary of sgRNA id and count
    {id:0, id:0, ...}
    unmapDic:
    A dictionary of unmapped read and count
    {}
    sgRNADic
    {sequence:id} dictionary
    datastat
    Statistics of datasets ({key:value})
    prefix
    prefix seq upstream variable region
    suffix
    suffix seq downstream variable region
    Return value
    -----------
    datastat
    a dictionary structure of statistics
    mapDic:
    A dictionary of sgRNA id and count
    {id:count, id:count, ...}
    unmapDic:
    A dictionary of unmapped read and count
    {seq:count, seq:count, ...}
    '''
    nline=0
    logging.info('Parsing file '+filename+'...')
    nreadcount=0
    if filename.upper().endswith('.GZ'):
        import gzip
        openobj=gzip.open(filename,'rt')
    else:
        openobj=open(filename)
    nsynerror=0
    nunknown=0
    for line in openobj:
        # line=line.encode('utf-8')
        nline=nline+1
        if nline%1000000==1:
            logging.info('Processing '+str(round(nline/1000000))+ 'M lines..')
        if nline%4 == 2:
            nreadcount+=1
            read=line.strip()
            sgRNAid=check_read(read,prefix,suffix,args.variable_region_len,sgRNADic)
            if sgRNAid!='unknown' and sgRNAid!='synthesis error':
                mapDic[sgRNAid]=mapDic[sgRNAid]+1
            elif sgRNAid=='synthesis error':
                nsynerror+=1
                # save unmapped file
                if args.unmapped_to_file:
                    if read not in unmapDic:
                        unmapDic[read]=0
                    unmapDic[read]=unmapDic[read]+1
            elif sgRNAid=='unknown':
                nunknown+=1
                # save unmapped file
                if args.unmapped_to_file:
                    if read not in unmapDic:
                        unmapDic[read]=0
                    unmapDic[read]=unmapDic[read]+1
            else:
                logging.info('incorrect read category %d'%(nline))
    openobj.close()
    print ('%s pretreatment finalized!'%(filename))
    # calculate statistics
    datastat['reads']=nreadcount
    nmapped=0
    nrdcnt=[]
    nzerosgRNA=0
    for (k,v) in mapDic.items():
        if v>0:
            nmapped+=v
            nrdcnt+=[math.log(v+1.0)]
        elif v==0:
            nzerosgRNA+=1
            nrdcnt+=[math.log(0.0+1.0)]
    logging.info('mapped:'+str(nmapped))
    datastat['mappedreads']=nmapped
    datastat['unmap due to syn error']=nsynerror
    datastat['unmap due to unknown source']=nunknown
    datastat['totalsgRNAs']=len(sgRNADic)
    datastat['zerosgRNAs']=nzerosgRNA
    datastat['giniindex']=CRISPRsgRNA_gini(nrdcnt)
    #return ctab
    return (mapDic,unmapDic)


# print partion!

def CRSIPRsgRNA_printdict(dict0,unmapdict0,args,ofile,ounmappedfile,datastat,sep='\t'):
    '''
    Write the table count to file
    '''
    allfastq=args.fastq.split(',')
    nsample=len(allfastq)
    slabel=[datastat[f]['label'] for f in allfastq]
    # print header
    print('sgRNA'+sep+sep.join(slabel),file=ofile)
    print('unmapped read'+sep+sep.join(slabel),file=ounmappedfile)
    # print mapped sgRNAs
    for (k,v) in dict0.items():
        print(k+sep+sep.join([str(x) for x in v]),file=ofile)
    # print unmapped reads
    for (k,v) in unmapdict0.items():
        print(k+sep+sep.join([str(x) for x in v]),file=ounmappedfile)

def CRSIPRsgRNAcount_printdict(dict0,args,ofile,datastat,sgRNA_geneDic,sep='\t'):
    '''
    Write the Normalize count table to file
    '''
    
    allfastq=args.fastq.split(',')
    nsample=len(allfastq)
    slabel=[datastat[f]['label'] for f in allfastq]
    # print header
    print('sgRNA'+sep+'Gene'+sep+sep.join(slabel),file=ofile)
    for (k,v) in dict0.items():
        print(k+sep+sgRNA_geneDic[k][0]+sep+sep.join([str(x) for x in v]),file=ofile)



def CRISPRsgRNA_printstat(args,datastat):
    '''
    Write data statistics to PDF file
    '''
    Xaxis=[]
    Yaxis=[]
    for (k,v) in datastat.items():
        logging.info('Summary of file '+k+':')
        for (v1,v2) in v.items():
            logging.info(str(v1)+'\t'+str(v2))
    # write to table
    crv=VisualRCount()
    crv.setPrefix(args.output_prefix)
    for (fq, fqstat) in datastat.items():
        crv.fastqfile+=[fq]
        if 'label' in fqstat:
            crv.fastqlabels+=[fqstat['label']]
            Xaxis.append(fqstat['label'].split('/')[-1])
        else:
            crv.fastqlabels+=['NA']
        if 'reads' in fqstat:
            crv.reads+=[fqstat['reads']]
        else:
            crv.reads+=[0]
        if 'mappedreads' in fqstat:
            crv.mappedreads+=[fqstat['mappedreads']]
        else:
            crv.mappedreads+=[0]
        if 'unmap due to syn error' in fqstat:
            crv.unmapduetosynthesis+=[fqstat['unmap due to syn error']]
        else:
            crv.unmapduetosynthesis+=[0]
        if 'unmap due to unknown source' in fqstat:
            crv.unmapunknown+=[fqstat['unmap due to unknown source']]
        else:
            crv.unmapunknown+=[0]
        if 'zerosgRNAs' in fqstat:
            crv.zerocounts+=[fqstat['zerosgRNAs']]
        else:
            crv.zerocounts+=[0]
        if 'giniindex' in fqstat:
            crv.gini+=[fqstat['giniindex']]
            Yaxis.append(fqstat['giniindex'])
        else:
            crv.gini+=[0.0]
    n_groups=len(Xaxis)
    index = np.arange(n_groups)
    plt.bar(index, Yaxis, 0.2,alpha=0.3,color='r',label='Gini index')
    plt.ylim(0,1)
    plt.xlabel('Library')
    plt.ylabel('Gini score')
    plt.title('Gini Scores by Library')
    plt.yticks(size=6)
    plt.xticks(index , np.array(Xaxis),rotation=300,size=6)
    plt.legend()
    plt.subplots_adjust(left=0.18, wspace=0.25, hspace=0.25,bottom=0.23, top=0.9)
    plt.savefig('%s_Libray_Gini_Score.png'%(args.output_prefix),dpi=1000)
    plt.clf()
    #
    crv.writeCountSummary()
    # write to TXT file
    crv.writeCountSummaryToTxt(args.output_prefix+'.countsummary.txt')


def CRISPRsgRNA_main(args):
    """
    Main entry for CRISPRsgRNA count module
    """
    # check arguments
    CRISPRsgRNA_checkargs(args)
    # check the listed files
    # listfq=args.fastq.split(',')
    listfq=args.fastq.split(',')
    nsample=len(listfq)
    prefix_nucl=args.prefix_nucl
    suffix_nucl=args.suffix_nucl
    datastat={}
    # check labels
    alllabel=args.sample_label
    if alllabel=='':
        slabel=['sample'+str(x) for x in range(1,nsample+1)]
    else:
        slabel=alllabel.split(',')
    for i,fi in enumerate(listfq):
        datastat[fi]={}
        datastat[fi]['label']=slabel[i]
    # process sgRNA dicts
    sgRNADic={}
    sgRNA_geneDic={}
    if args.list_seq is not None:
        (sgRNA_geneDic,sgRNADic,sgRNALst)=CRISPRsgRNA_checklists(args)
    mapdict={}
    unmapdict={}
    # go through the fastq files
    for filename in listfq:
        dict0={sgRNAid:0 for sgRNAid in sgRNALst}
        undict0={}
        CRISPRsgRNA_processonefile(filename,args,dict0,undict0,sgRNADic,datastat[filename],prefix_nucl,suffix_nucl)
        CRISPRsgRNA_mergedict(mapdict,dict0)
        #CRISPRsgRNA_mergedict(unmapdict,undict0)
    # write to file
    # generate the files
    list_files=args.output_prefix.split('/')
    output_dir=''
    for i in range(len(list_files)-1):
        output_dir=output_dir+list_files[i]+'/'
    os.system('mkdir -p %s' %(output_dir))
    fl=args.output_prefix+'.count.txt'
    os.system('cd %s'%(output_dir))
    os.system('cat /dev/null > %s'%(fl))
    ofilel=open(fl,'w')
    if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
        ounmappedfilel=open(args.output_prefix+'.unmapped.txt','w')
    else:
        ounmappedfilel=None
    CRSIPRsgRNA_printdict(mapdict,unmapdict,args,ofilel,ounmappedfilel,datastat,sep='\t')
    ofilel.close()
    if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
        ounmappedfilel.close()
    os.system('rm %s.unmapped.txt'%(args.output_prefix))
    f2=args.output_prefix+'.normalizeCount.txt'
    os.system('cat /dev/null > %s'%(f2))
    ofile2=open(f2,'w')
    medalldict=normalizeCount(mapdict)
    CRSIPRsgRNAcount_printdict(medalldict,args,ofile2,datastat,sgRNA_geneDic,sep='\t')
    CRISPRsgRNA_printstat(args,datastat)
    return 0


if __name__ == '__main__':
    try:
        args=CRISPRsgRNA_parseargs()
        CRISPRsgRNA_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)
