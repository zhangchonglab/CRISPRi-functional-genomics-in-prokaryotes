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

5. Operon file. CRISPRi works at the transcription level, due to the unique structure of polycistronic operons in prokaryotic genomes, it is hard to figure our the true phenotype-associated genes when multiple genes in one polycistronic operon responses. To cope with this problem, we designed in this subpackage to reorganize the gene level statistics at the operon level as an option. If you are not interested in this step, please ignore this file and no need to prepare it. 
The file has header line and is consisted of three columns: operon id, operon name and the genes in the operon. Tab is used as delimiter, and the genes are separated by commas in each line. We recommand to organize genes in one polycistronic operon according to the order from upstream to downstream. Gene names should be consistent with those in sgRNA-library file and sgRNA-position file.
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
Step 3: Edit the configure file, which is used to set all the necessary parameters for the design. The default parameters are given in the original example configure file. Generally, lines beginning with '#' is the anootation; other lines start with one word seperated with the following by a tab delimiter. All files (subdirectories) mentioned in this configure file should be located under a common working directory with all .py scripts.

prefix: Prefix used for naming of all output files.

fastqpath: directory under which there are all .fastq files

fastq: the name of .fastq files to be processed under 'fastqpath'. The files can be fastq.gz or .fastq format. Multiple .fastq file names are seperated by ',' (comma).

forward_prefixseq: several (4-10) nucleotides upstream (usually the last several nucleotides in promoter) the variable region (protospacer) of sgRNA used to specify and cut the variable region from the sequencing read. These nucleotides should be within the prepared sequencing library.

forward_suffixseq: several (4-10) nucleotides downstream (in the Cas9 binding motif) the variable region (protospacer) of sgRNA used to specify and cut the variable region from the sequencing read. These nucleotides should be within the prepared sequencing library.

sample-label: the label for each sequencing library seperated by comma. The order of the label corresponds to the order of the .fastq file names specified by the 'fastq' option. Hence, the total number of label should be tha same as that of file names. See example configure file to understand this intuitively. Note that the labels specified here should be the same as the library name defined in the experiment design file (see above). 

sgrna-len: number of nucleotides of the variable (protospacer) region of the sgRNA. It is specified in the library design. Note that the length specified here should be consistent with that of the sgRNA-library file (see above). default: 20.

list-seq: the name of the sgRNA-library file (see above).

experiment_configure: the name of the experiment design file (see above).

control_setting: sgRNAs used as control to calculate the statistics of the gene-phenotype association. Two options: 'NC' and 'all'. In the case where negative control sgRNAs were included in the synthetic library during the selection experiment, 'NC' is recommanded. In other cases where no negative control sgRNAs are available, 'all' should be specified. We highy recommand to include negative control sgRNAs during the experiment and data analysis, because this option significantly improves the statistical robustness of the result. Note that when 'NC' is specified here, negative control sgRNAs should be included in the sgRNA-library file (see above).

Step 4：Run the script

Examples
--------
python CRISPRscreen_main.py configure.txt
