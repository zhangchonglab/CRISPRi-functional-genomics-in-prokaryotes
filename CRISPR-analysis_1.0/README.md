# synthetic sgRNA library based functional genomics for prokaryotes: NGS data processing part

## What is this?
This python script collection is one of the two the software subpackages of CRISPRi functional genomics method for the prokaryotes, used for the NGS data processing and result visualization. The basic description of this program can be found at BioRxiv: https://doi.org/10.1101/129668. Please cite this paper or subsequent peer-reviewed publication if this program is useful to your work.

This script collection is user-friendly for experimental microbiologists with no or limited programming expertise. Generally, the user only need to download the script, edit a configure file to set several parameters needed for sgRNA design, and type in one command line in a Linux environment to initiate the design process. The output includes the statistics about result at sgRNA, gene and operon level of each stressed condition. Meanwhile, the Gini index profile of each library, the agreement of biological replicates, fitness distribution of negative sgRNA, FPR-score curve and diverse visualization files are also presented.

## General description of the algorithm and experiment
The synthetic sgRNA plasmid library designed by the sgRNA-design subpackage can be transformed into prokaryotic cells expressing dCas9 protein. The resulting cell library can be subjected to stressed and control condition, and the plasmids after selection can be extracted and prepared for NGS library. With this protocol, we have a series of NGS data (.fastq) for each stressed and control condition. Each one can have one or more biological replicates. About the detailed protocol of the method, see https://doi.org/10.1101/129668 in BioRxiv. This program is used to convert these .fastq data to the gene-phenotype association profile.

## How to use it
### Step 1：Installation
1. Install Python version 2.7 or above 
2. Install Scipy version 0.19.1 or above
3. Install Matplotlib version 2.0.2 or above
4. Install Numpy version 1.13.1 or above

### Step 2：Prepare the necessary files.
All these files (or subdirectories) should be organized under a common working directory together with the all .py scripts. The structure of the working directory is shown as below as the checklist:
under a common working_dir/[example illustration of files under the working directory](./image/files_prepared_before_data_processing.png)
Please check the example files post at GitHub, which are described as below.

#### File 1: NGS files (.fastq or .fq extension) under one directory (example_data/all files)
Note: Please try to keep the name of each file meaningful but as simple as possible, thus simplifying the preparation of configure file as described below.

#### File 2: sgRNA library file (example_library.csv)
The sgRNA library file is at .csv formate **containing the header line**, in which three are three columns that are in order of id, sequence and gene respectively. **Use comma as delimiter**. 
If negative control sgRNAs are within this synthetic library, name them NCx and assign '0' at 'gene' column of these sgRNAs. This file can be found as an output of the library design subpackage. It should be noted that **Note that -, _ and ' '(space) should be eliminated from any id name. Avoid id like 'super-sgRNA', 'super_sgRNA' or 'super sgRNA'**.

id|sequence|gene
--|--------|----
sgRNA1|ATCCCCCCCCCCGGGGG|recA
sgRNA2|ACTGCCCCGGGGCCCCC|recA
NC1|ACACACACACACACACACAC|0
NC2|TGTGTGTGTGTGTGTGTGTG|0
...|...|...

#### File 3: sgRNA position file (example_coding_region_position.txt)
Flat file of sgRNA position (relative location of sgRNA in the coding region) information in gene **without header line using tab as delimiter**.
The file is without header line containing three columns that are in order of gene name, sgRNAid and the relative position of sgRNA in the gene. The name of sgRNAid and gene should be compatible with the sgRNA library file (File 2). Actually, you can also find this file as output of the library design subpackage.

rsmE|rsmE_9|0.0122950819672
----|------|---------------
rsmE|rsmE_10|0.0136612021858
rsmE|rsmE_11|0.0150273224044
rsmE|rsmE_12|0.016393442623
rsmE|rsmE_25|0.0341530054645
acnA|acnA_384|0.143497757848
acnA|acnA_395|0.147608370703
acnA|acnA_441|0.164798206278
acnA|acnA_459|0.171524663677
acnA|acnA_477|0.178251121076
...|...|...

#### File 4: operon file (example_operon.txt) **This is optional**
CRISPRi works at the transcription level, due to the unique structure of polycistronic operons in prokaryotic genomes. It is hard to figure out the true phenotype-associated genes coping with multiple genes in one polycistronic operon. To address this problem, during the design of this package, we reorganize the gene level statistics at the operon level as an option. If you are not interested in this step or in other cases that your microorganism do not have a available operon file, please ignore this and no need to prepare it. 
**The file has a header line and uses tab as delimiter.** It is consisted of three columns: operon id, operon name and the genes in the operon. **Genes in one polycistronic operon are separated by comma**. If one gene is located at multiple operons, it is ok to just list all of them. We recommend to organize genes in one polycistronic operon according to the order from upstream to downstream. Gene names should be consistent with those in sgRNA library file (File 2) and sgRNA position file (File 3).

operonid|operon	genes
|-------|------------
KO04087|dapE,ypfN
KO04086|yffB,dapE,ypfN
KO04089|bcp,gcvR
KO04735|ykfH,ykfF,yafX,ykfI,ykfG,yafW
KO04731|dinJ,yafQ
KO04956|ssuB,ssuE,ssuA,ssuD,ssuC
KO05736|damX,aroB,gph,trpS,dam,rpe,aroK
...|...

#### File 5: experiment design file (example_experiment_configure.txt)
For each phenotype to be studied, we need one selective and one control condition, respectively. This file is used to define the role of each screening experiment (NGS data). **This file has one header line and uses tab as delimiter. The header line should be organized in a format as below:**

Library/Condition|initial|stress1|control1|...|stressX|controlX|...
-----------------|-------|-------|--------|---|-------|--------|---

**Each row refers to a NGS library and each column refers to a condition. '1' indicates the association between library and condition, whereas the program will skip the item set as ‘0’.** All libraries under one common condition (in one column) are regarded as biological replicates and read count for each sgRNA of these libraries are averaged as geometric mean. One library can be associated with multiple conditions. **At least and only one library should be assigned as initial condition.** Initial library will be used to exclude sgRNAs with poor representation from further analysis to ensure statistical robustness based on customized threshold (see below in configure file). Usually, use the library before selection as the initial library. Initial library can be also used as control condition for a particular phenotype to be studied.

Library/Condition|initial|stress1|control1
-----------------|-------|-------|--------
dCas9R1|0|1|0
dCas9R2|0|1|0
NCR1|0|0|1
NCR2|0|0|1
slib1-10mixNC|1|0|0

#### File 6: naming file (example_naming_configure.txt)
To name output files related to different studied phenotypes, we design this naming file to give each phenotype a name.
**This file has no header line and uses tab as delimiter.***
The file has two columns, while the first is a key pointing to File 5 (the experiment design file) and the second defines the name for each phenotype. Hence, the items in the first column should be like 'stressX', thus compatible with File 5.

stress1|essential
-------|---------
stress2|whatever you like
stress3|whatever you like
...|...

///////////////////////////////////////////////////////////////
### Step 3: Set up the configure file
Edit the configure file, which is used to set all the necessary parameters and path towards necessary files. The default parameters are given in the original example configure file. Generally, lines beginning with '#' is the annotation; other lines start with one word (name of this parameter) separated with the following (setting of this parameter) by a tab delimiter.
 
prefix: Prefix used for naming of all output files, keep it simple without any ‘_’ or ‘-’. For example, ‘myscreen2017’ is fine.

fastqpath: directory under which all .fastq files are located. See step 2 and part 1 above.

fastq: the name of .fastq files to be processed under 'fastqpath'. The files can be fastq.gz or .fastq format. Multiple .fastq file names are separated by ',' (comma).

forward_prefixseq: several (4-10) nucleotides upstream (usually the last several nucleotides in promoter) the variable region (protospacer) of sgRNA used to specify and cut the variable region from the sequencing read. These nucleotides should be located in the range of the PCR product during sequencing library preparation.

forward_suffixseq: several (4-10) nucleotides downstream (in the Cas9 binding motif) the variable region (protospacer) of sgRNA used to specify and cut the variable region from the sequencing read. These nucleotides should be located in the range of the PCR product during sequencing library preparation.

sample-label: the label for each sequencing library separated by comma. The order of the label should correspond to the order of the .fastq file names specified by the 'fastq' option. Hence, the total number of label should be the same as that of file names. See example configure file to understand this intuitively. Note that the labels specified here should be the same as the library name defined in the experiment design file (see above, Step 2, Part 5).

sgrna-len: number of nucleotides of the variable region (protospacer, also the nucleotides between prefixseq and suffixseq) of the sgRNA. It is specified in the library design. Note that the length specified here should be consistent with that of the sgRNA-library file (see above). default: 20.

list-seq: the name of the sgRNA-library file (see above, Step 2, Part 2).

experiment_configure: the name of the experiment design file (see above, Step 2, Part 5).

name_configure: the name of the naming file for each stressed condition (see above, Step 2, Part 6)

control_setting: sgRNAs used as control to calculate the statistics of the gene-phenotype association. Two options: 'NC' or 'all'. In the case where negative control sgRNAs (sgRNA targeting nowhere in the genome, set as gene=‘0’, see Step 2, Part 2) were included in the synthetic library during the selection experiment, 'NC' is recommended. In other cases where no negative control sgRNAs are available, 'all' should be specified. We highly recommend to include negative control sgRNAs during the experiment and data analysis, because this option significantly improves the statistical robustness of the experiment. Note that when 'NC' is specified here, negative control sgRNAs should be included in the sgRNA-library file (see Step 2, Part 2).

FDR_threshold: FDR (False discovery rate) threshold to call significant phenotype-associated genes. Default: 0.01. For details to calculate FDR for each gene-phenotype association, see the paper.

hit_gene_calling: method to call hit gene associated with particular phenotype. Two options: 'position' or 'all'. For ‘position’ option, the program search sgRNA fitness belonging to a gene one by one according to their position within the gene ORF (start from those proximal to start codon) (specified by sgRNA-position.txt, Step 2, Part 3), and used the sgRNA subset with smallest FDR value. For ‘all’ option, all sgRNAs belonging to the gene were used. Due to the fact we found that sgRNAs targeting to the 5' of ORF exhibited better repression activity, 'position' method is more recommended. For details, see the paper.

gene_sgRNA_position: name of the flat file of sgRNA position (see above, Step 2, Part 3).

Operon_gene_List: name of the operon file (see above, Step 2, Part 4). This one is optional, if you do not need it, just leave it blank.

After Step 2 and 3, check your working directory. It should looks like below:
under the working_dir/
*.py rawdata_dir/*.fq library.csv sgRNA-position.txt operon.csv(optional) experiment.txt naming.txt configure.txt

=============================================
### Step 4：Run the script

Open the command line window, cd to the working directory and run the pipeline.
--------
cd working_directory (the directory containing all necessary files mentioned above and .py scripts)
python CRISPRscreen_main.py configure.txt

If you want to test your environment, we have a toy example for testing purpose. cd to the working_dir, type in: 
python CRISPRscreen_main.py configure.txt

///////////////////////////////////////////////////////////////
## Output files
The output files would be all organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory.
