# synthetic sgRNA library based functional genomics for prokaryotes: sgRNA library design package

## What is this?
This python script collection is one of the two the software subpackages of CRISPRi functional genomics method for the prokaryotes, used for genome-wide or focused sgRNA library design. The basic description of this program can be found at BioRxiv: https://doi.org/10.1101/129668. Please cite this paper or subsequent peer-reviewed publication if this program is useful to your work.

This script collection is user-friendly for experimental microbiologists with no or limited programming expertise. Generally, the user only need to download the script and several standard files (genome, annotation, etc), edit a configure file to set several parameters needed for sgRNA design, and type in one command line in a Linux environment (for example, terminal in MacOS) to initiate the design process. The output includes the .fasta file for sgRNA and statistics-describing tables as well as figures to intuitively overview the design result.

Seqmap, developed by Jiang, H., Wong, W.H. (2008) Bioinformatics, 24(20) is used in our program for off-target identification and elimination. We have included an executable file of seqmap in the delivered zip file. If you want to use the updated version of seqmap, please follow the below instructions:
1. Download the Seqmap source file (for all platform) from http://www-personal.umich.edu/~jianghui/seqmap/ and unzip it.
2. Compile the source code following the README file under the unzipped directory of seqmap. 
3. Copy and paste the resulted seqmap executable file to the working directory.

## How to use it?
### Step 1: Installation
1. Install Python version 2.7 or above
2. Install Matplotlib version 2.0.2 or above
3. Install Numpy version 1.13.1 or above

### Step 2: Download or construct the necessary files. 
This package provides two options for library design: genome-wide or focused sgRNA library. 
1. genome-wide library
For one microorganisms, three (or five if RNA-coding genes are also of interest) standard files are needed to design the genome-scale sgRNA library for all protein-coding genes by the genome. All these files are accessed at the FTP site of NCBI: ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/
For example, for *E. coli* K12 MG1655 strain, open the hyperlink above, find the directory for this strain:
ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/ [like here](./image/NCBI_refseq_E.coli.png)

2. focused library
The detailed description about these files can be found at Step 3 (targetFasta, indexFile, genome).
Please download the above-mentioned files from the ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/. You can find the relevant files under the directory of specific organism of your own interest. Please download or construct these files into the working directory where the scripts and configure file are located.

### Step 3: Edit the configure file, which is used to set all the necessary parameters for the design. 
The default parameters are given in the original configure file under the unzipped directory. Generally, the program processes the gene one by one, search for all the possible sgRNAs for this gene meeting the requirement given in the configure file, eliminate the sgRNAs with potential off-target site(s), and store all the qualified sgRNAs. The brief description of all the parameters is given below. For more details, please check our paper (BioRxiv: https://doi.org/10.1101/129668).

**ORFcutoff**: The sgRNA location within the gene coding region (ORF 5’=0.0, ORF 3’=1.0) (default=0.05, real number belonging to (0.0,1.0)). We found that the sgRNAs exhibited higher activities locating with the first 5% of ORF. Hence, by default, the sgRNA is designed at this region as many as possible.

sgRNA_number: The number of sgRNA you want to design for each gene (default=10, positive integer accepted). According to the description of our paper, we found that in E. coli, 5-10 sgRNAs per gene is enough for pooled screen based functional genomics analysis. However, more sgRNAs provide benefit (stronger statistical significance) for genes with only moderate phenotypic effect.

GCcontentMin: The minimal GC content of spacer region (percentage) (default=20, positive integer between 0~100 is accepted).

GCcontentMax: The maximum GC content of spacer region (percentage) (default=80, positive integer between 0~100 and > GCcontentMin is accepted). 

In previous reports about dCas9 based CRISPRi system, GC content of spacer region is correlated with sgRNA activity only in extreme GC content (reduce the activity). Hence, we suggested the above threshold. In situations of genome with relative low or high GC content, we suggested to adjust the threshold to (10,90). 

off_threshold: Set the off target penalty threshold (default=20), sgRNAs with potential off-target site carrying penalty score lower than the threshold will be eliminated. For the detailed description of the scoring method, please check the paper. Briefly, we suggested off_threshold >= 20 for library design. In extreme situation where more sgRNAs are desired, the threshold can be decreased to 10.

strand: Whether the sgRNA is designed targeting (binding) to the template or nontemplate stand of a coding gene (default=nontemplate, nontemplate or template is accepted). It is suggested by previous reports that dCas9 based CRISPRi system used in this work exhibited higher activity when targeting to non-template strand in ORF region.

negative: Choose whether to design negative control sgRNAs (sgRNA with no significant target across the genome, which is used as negative control in pooled screen to determine the phenotypic effect and statistical power of each gene) for the experiment(yes or no， default=yes). We highly recommend to include the negative control sgRNAs in the pooled screen. For the description of negative control sgRNA usage, please see our paper.

negative_number: The number of negative control sgRNA you want to design for the experiment.If negative option is no, select 0 for this option. The default is 400.

targetFasta: Gene sequence file for genes of interest in .fasta format. It can be downloaded from NCBI ftp site for all protein-coding gene (.ffn suffix) or all RNA-coding genes (.frn suffix) encoded by a genome of particular microorganisms. The customized .fasta file for a focused library design can be extracted from the abovementioned files. It is noted that the naming of each gene in the customized .fasta file should follow that of .ffn or .frn file from NCBI.

indexFile: The gene sequence annotation file. (.ptt suffix for protein-coding genes and .rnt suffix for RNA-coding genes).

genome: The genome file (.fna suffix) used for off-target check. 

prefix: Prefix used for naming of all output files.

### Step 4: After editing the configure file as introduced above, cd to the working directory where the scripts and the necessary files (check list: configure file, files for targetFasta, indexFile and genome, all python scripts, seqmap executable file) are located. Type in the command line below:

python sgRNA_desgin_main.py configure.txt

## Output description
All the output files will be located in a directory named after the prefix given in configure.txt, under which the files about negative control sgRNA are under a subfolder named after negative.
The prefix.fasta file contains the the target gene sequences.
The prefix.txt is the .txt format of output.fasta file.
The prefix.N20.fasta contains the sgRNA library sequences (N20).
The prefix.N20.txt is the .txt format of output.N20.fasta file.
The prefix.N20NGG.fasta contains the sgRNA library sequence with PAM in the relevant genome.
The prefix.N20NGG.txt is the .txt format of output.N20NGG.fasta file.
The prefix.sgRNA_statistics.txt file contains the position information within relevant ORF and the GC content of each designed sgRNA.
The prefix.gene_statistics.txt file contains the length and the designed sgRNA numbers of each gene.
Two .png figures use histogram to summarize the basic information of gene (sgRNA number vs. gene) and sgRNA (position vs. sgRNA) statistics.

Under the negative subdirectory:
The prefix_N20_NC_passed.fasta contains the negative control sgRNA sequences.
The prefix_N20_NC_passed.txt is the .txt format of the output_N20_NC_passed.fasta.

We also included example files of Corynebacterium glutamicum ATCC 13032 (RefSeq NC_003450). After unzipping the download, please directly come to Step 4 to run an example process as an initial test (79 C. glutamicum genes, around 30 min running).

A typical genome-wide design process needs around 24h running on a typical laptop (for example, 2.6 GHz processor, 8 GB memory). Please keep the power on until the process ends automatically. 
