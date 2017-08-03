======================
synthetic sgRNA library based functional genomics for prokaryotes: NGS data processing part
======================

# ///////////////////////////////////////////////////////////////
What is this?
This python script collection is one of the two the software subpackages of CRISPRi functional genomics method for the prokaryotes, used for the NGS data processing and result visualization. The basic description of this program can be found at BioRxiv: https://doi.org/10.1101/129668. Please cite this paper or subsequent peer-reviewed publication if this program is useful to your work.

This script collection is user-friendly for experimental microbiologists with no or limited programming expertise. Generally, the user only need to download the script, edit a configure file to set several parameters needed for sgRNA design, and type in one command line in a Linux environment to initiate the design process. The output includes the statistics about result at sgRNA, gene and operon level of each stressed condition. Meanwhile, the Gini index profile of each library, the agreement of biological replicates, dilution distribution of negative sgRNA, FDR-p value curve raw data and diverse visualization files are also presented.

# ///////////////////////////////////////////////////////////////
General description of the algorithm and experiment
The synthetic sgRNA plasmid library designed by the sgRNA-design subpackage can be transformed into prokaryotic cells expressing dCas9 protein. The resulting cell library can be subjected to stressed and control condition, and the plasmids after selection can be extracted and prepared for NGS library. With this protocol, we have a series of NGS data (.fastq) for each stressed and control condition. Each one can have one or more biological replicates. About the detialed protocol of the method, see https://doi.org/10.1101/129668 in BioRxiv. This program is used to convert these .fastq data to the gene-phenotype (stree) association profile.


Step 1：Installation
============
1. Install Python version 2.7 or above 
2. Install BioPython version 1.57 or above
3. Install Scipy version 0.17.1 or above
4. Install Matplotlib version 2.0.0 or above
5. Install Numpy version 2.0.0 or above

Step 2：Prepare the necessary files.
All these files (subdirectories) should be organized under a common working direcotry together with the all .py scripts.
You can find all example files corresponding to the belowmentioned files.
1. All .fastq files under one directory.

2. The sgRNA-library file that is a csv formate file containing the header line, in which three are three columns that are in order of id, sequence and gene respectively, and each line uses a comma as delimiter. 
If negative control sgRNAs are within this synthetic library, name them whatever you like (NCx for example) and assign '0' at 'gene' column of these sgRNAs. 
(For example :
	id,sequence,gene
	sgRNA1,ATCCCCCCCCCCGGGGG,recA
	sgRNA2,ACTGCCCCGGGGCCCCC,recA
	NC1,AAAAAAAAAAAAAAAAAAAA,0)
This file can be found as an output of the library design subpackage. It should be noted that '-' should be eliminated from any id.

3. The experiment design file. This file is used to distinguish between the initial (before selection), stressed and control (after selection) conditions with tab as delimiter. 
Each row refers to a seqeuncing library and each column refers to a contidtion. '1' indicates the association between library and condition. All libraries under one common condition are regarded as biological replicates and read count for one sgRNA of these libraries are averaged as geometric mean. One library can be associated with multiple conditions. At least one library should be assigned as initial condition (before selection). Initial library will be used to exclude sgRNAs with poor representation from further analysis to ensure statistical robustness based on customized threshold (see below). Usually, keep cell library before selection to prepare initial seqeuncing library. Otherwise, you can use particular library under control condition as initial library (assume no sgRNA inhibits or improves cell growth together with dCas9 in this control condition).   
(For example:
	Library/Condition	initial	stress1	control1	stress2	control2
	M_LB_C1_R1_1	0	0	1	0	0
	M_LB_C1_R2_1	0	0	1	0	0
	M_LB_C2_R1_1	0	0	0	0	1
	M_LB_C2_R2_1	0	0	0	0	1
	Minimal_Bf_1	1	0	0	0	0	
	MOPS_C1_R1_1	0	1	0	0	0
	MOPS_C1_R2_1	0	1	0	0	0
	MOPS_C2_R1_1	0	0	0	1	0
	MOPS_C2_R2_1	0	0	0	1	0)
	
4. Flat file of sgRNA position (relative location of sgRNA in ORF of relevant target gene, see paper for details) information in gene, this file is in current version for hit-gene calling because generally sgRNAs locating at 5' of ORF exhbited better knockdown activity. The file is without header line containing three columns that are in order of gene name, sgRNAid and the relative position of sgRNA in the gene，and each line uses a tab as a delimiter.
(For example :
	rsmE	rsmE_9	0.0122950819672
	rsmE	rsmE_10	0.0136612021858
	rsmE	rsmE_11	0.0150273224044
	rsmE	rsmE_12	0.016393442623
	rsmE	rsmE_25	0.0341530054645
	acnA	acnA_384	0.143497757848
	acnA	acnA_395	0.147608370703
	acnA	acnA_441	0.164798206278
	acnA	acnA_459	0.171524663677
	acnA	acnA_477	0.178251121076)
This file can be found as an output of the library design subpackage.

5. Operon list file.This file is used to statistic the sgRNA on the level of Operon.If you do not want to do it ,you do not prepare this file. The file with title is made consist of three columns, and they are Operon numbering, Operon name and the genes that are contained by Operon. The file formate
 is that each line uses a tab as a delimiter, and the genes are separated by commas in each line.
(For example:
		koid	name	op	
		KO04087	TU-1463	dapE,ypfN
		KO04086	TU-1462	yffB,dapE,ypfN
		KO04089	TU-1466	bcp,gcvR
		KO04735	TU-2623	ykfH,ykfF,yafX,ykfI,ykfG,yafW
		KO04731	TU-2616	dinJ,yafQ
		KO04956	TU-2970	ssuB,ssuE,ssuA,ssuD,ssuC
		KO05736	TU-4212	damX,aroB,gph,trpS,dam,rpe,aroK
		KO04558	TU-2303	ubiA,ubiC
		KO04884	TU-2858	kdpF,kdpA,kdpC,kdpB
		KO04554	TU-2299	yjbF,yjbE,yjbH,yjbG
		KO06171	TU-Exp-0624	sfsA,dksA)

Step 3：Set up the configure file
# ///////////////////////////////////////////////////////////////
Step 3: Edit the configure file, which is used to set all the necessary parameters for the design. The default parameters are given in the original configure file under the unzipped directory. Generally, the program processes the gene one by one, search for all the possible sgRNAs for this gene meeting the requirement given in the configure file, eliminate the sgRNAs with potential off-target site(s), and store all the qualified sgRNAs. The brief description of all the parameters is given below. For more details, please check our paper (BioRxiv: https://doi.org/10.1101/129668).

ORFcutoff: The sgRNA location within the gene coding region (ORF 5’=0.0, ORF 3’=1.0) (default=0.05, real number belonging to (0.0,1.0)). We found that the sgRNAs exhibited higher activities locating with the first 5% of ORF. Hence, by default, the sgRNA is designed at this region as many as possible.

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

Step 4：Run the script

Examples
--------
python CRISPRscreen_main.py configure.txt
