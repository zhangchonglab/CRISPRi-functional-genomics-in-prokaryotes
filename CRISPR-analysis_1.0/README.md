# synthetic sgRNA library based functional genomics for prokaryotes: NGS data processing package

## What is this?
This python script collection is one of the two the software subpackages of CRISPRi functional genomics method for the prokaryotes, used for the NGS data processing and result visualization. The basic description of this program can be found at BioRxiv: https://doi.org/10.1101/129668. Please cite this paper or subsequent peer-reviewed publication if this program is useful to your work.

This script collection is user-friendly for experimental microbiologists with no or limited programming skills. Generally, the user only need to download the script, edit a configure file to set several parameters needed for sgRNA design, and type in one command line in a Linux environment to initiate the process. The output includes the statistics at sgRNA, gene and operon level of each studied phenotype. Meanwhile, diverse visualization files are also presented.

## General description of the algorithm and experiment
The synthetic sgRNA plasmid library designed by the sgRNA-design subpackage can be transformed into prokaryotic cells expressing dCas9 protein. The resulting cell library can be subjected to stressed and control condition to study the reponse of genes to phenotypes, and the plasmids after selection can be extracted and sequenced by NGS. About the detailed protocol of the method, see our paper. This program is used to convert NGS raw data to the gene-phenotype association profile. For the NGS raw data mapping step of this package, we referred to the algorithm and source code described by Li, W. et al. Genome Biol. 16, 281 (2015) with some modifications.

## How to use it?
### Step 1：Installation
1. Install Python version 2.7 or above
2. Install Scipy version 0.19.1 or above
3. Install Matplotlib version 2.0.2 or above
4. Install Numpy version 1.13.1 or above

### Step 2：Prepare the necessary files.
All these files (or subdirectories) should be organized under a common working directory together with the all .py scripts.
Please check the example files post at GitHub, which are described as below.

#### File 1: NGS files (.fastq or .fq extension) under one directory (see example_data/all files)
Note: Please try to keep the name of each file meaningful but as simple as possible. The files can also be compressed as .gz format.

#### File 2: sgRNA library file (see example_library.csv)
The sgRNA library file is at .csv formate **containing one header line**, in which there are three columns in order of id, sequence and gene respectively. **Use comma as delimiter**.
If negative control (NC) sgRNAs are within this synthetic library, name them NCx and assign '0' at 'gene' column of these sgRNAs. This file can be found as an output of our library design subpackage. It should be noted that **-, _ and ' '(space) should be eliminated from any id name. Avoid id like 'super-sgRNA', 'super_sgRNA' or 'super sgRNA'**.

sgRNAID|sgRNAseq|gene
--|--------|----
sgRNA1|ATCCCCCCCCCCGGGGG|recA
NC1|TGTGTGTGTGTGTGTGTGTG|0
...|...|...

#### File 3: sgRNA position file (see example_coding_region_position.txt)
Flat file of sgRNA position (relative location of sgRNA in the coding region) information in gene **without header line using tab as delimiter**.
The file contains three columns in order of gene name, sgRNAid and the relative position of sgRNA in the gene. The name of sgRNAid and gene should be compatible with the sgRNA library file (File 2). Actually, you can also find this file as output of our library design subpackage.

rsmE|rsmE_9|0.012
----|------|-----
rsmE|rsmE_10|0.014
0|NC1|0
...|...|...

#### File 4: operon file (see example_operon.txt) (optional)
During the design of this package, we reorganize the gene level statistics at the operon level as an option. If you are not interested in this step or if your microorganism does not have a available operon file, please ignore this. 
**The file has a header line and uses tab as delimiter.** It is consisted of two columns: operon id and the genes in the operon. **Genes in one polycistronic operon are separated by comma**. If one gene is located at multiple operons, it is ok to just list all of them. We recommend to organize genes in one polycistronic operon according to the order from upstream to downstream. Gene names should be consistent with sgRNA library file (File 2) and sgRNA position file (File 3).

operonid|operon	genes
|-------|------------
KO04087|dapE,ypfN
KO04089|bcp
...|...

#### File 5: experiment design file (see example_experiment_configure.txt)
For each phenotype to be studied, we need one selective and one control condition, respectively. This file is used to define the role of each NGS raw data. **This file has one header line and uses tab as delimiter. The header line should be organized in a format as below:**

Library/Condition|initial|stress1|control1|...|stressX|controlX|...
-----------------|-------|-------|--------|---|-------|--------|---

**Each stress and control pair defines a phenotype to be studied.**

**Each row refers to a NGS library and each column refers to a condition. '1' indicates the association between library and condition, whereas the program will skip the item set as ‘0’.** All libraries under one common condition (in one column) are regarded as biological replicates and read count for each sgRNA of these libraries are averaged as geometric mean. One library can be associated with multiple conditions. **At least and only one library should be assigned as initial condition.** Initial library will be used to exclude sgRNAs with poor representation from further analysis based on customized threshold (see below in configure file). Usually, use the library before selection as the initial library.

Library/Condition|initial|stress1|control1
-----------------|-------|-------|--------
dCas9R1|0|1|0
dCas9R2|0|1|0
NCR1|0|0|1
NCR2|0|0|1
plasmid|1|0|0

#### File 6: naming file (see example_naming_configure.txt)
To name output files related to different studied phenotypes, we design this naming file to give each phenotype a name.
**This file has no header line and uses tab as delimiter.***
The file has two columns, while the first is pointer to File 5 (the experiment design file) and the second defines the name for each phenotype. Hence, the items in the first column should be like 'stressX', thus compatible with File 5.

stress1|essential
-------|---------
stress2|whatever you like
stress3|whatever you like
...|...

### Step 3: Set up the configure file (see example_configure.txt)
The configure file is used to set all the necessary parameters and tell the program where to find necessary files. 
**This file is in a two-column format using tab as delimiter.** Each line starts with one word (name of one parameter) separated with the following (setting of this parameter) by a tab delimiter. We describe each parameter as below.
 
**prefix**: prefix used for naming of all output files, keep it simple without any ‘-’, ‘_’ and ‘ ’. For example, ‘screen20171001’ is fine.

**fastqpath**: directory under which all NGS raw data (.fastq or .fq extension) files are located. See File 1 above.

**fastq**: the names of NGS raw data (.fastq or .fq extension) files to be processed under 'fastqpath'. Multiple file names are separated by ',' (comma).

**forward_prefixseq**: several (4-10) upstream nucleotides (promoter) flanking the variable region (protospacer) of sgRNA used to specify and cut the variable region from the sequencing read. These nucleotides should be located in the the PCR product of NGS library.

**forward_suffixseq**: several (4-10) downstream nucleotides (Cas9 binding motif) flanking the variable region (protospacer) of sgRNA used to specify and cut the variable region from the sequencing read. These nucleotides should be located in the the PCR product of NGS library.

**sample-label**: the label for each NGS raw data file. The order of the label should corresponds to the order of the file names specified by the 'fastq' parameter. For simplicity, it is fine to use raw data file name without extension as label. Note that the labels specified here should be the same as the library name defined in the experiment design file (File 5).

**sgrna-len**: number of nucleotides of the variable region (protospacer, also the nucleotides between forward_prefixseq and forward_suffixseq) of the sgRNA. It is determined in the library design. Hence, the length specified here should be consistent with that of the sgRNA-library file (see above). default: 20.

**list-seq**: the name of the sgRNA library file (see above, Step 2, File 2).

**experiment_configure**: the name of the experiment design file (see above, Step 2, File 5).

**name_configure**: the name of the naming file (see above, Step 2, File 6)

**control_setting**: sgRNAs used as control to determine the experimental noise and calculate the statistics of the gene-phenotype association. **Two options: 'NC' or 'all'.** In the case where NC sgRNAs (sgRNA targeting nowhere in the genome, set as gene=‘0’, see Step 2, File 2) are included in the sgRNA library, 'NC' is recommended. In other cases, 'all' should be specified. **We highly recommend to include NC sgRNAs during the experiment and data analysis**, because this option significantly improves the statistical robustness of the experiment. Note that when 'NC' is specified here, NC sgRNAs should be included in the sgRNA-library file (see Step 2, File 2).

**FDR_threshold**: FDR (False discovery rate) threshold for hit gene calling. Default: 0.05.

**ReadsThreshold**: threshold of read count used to remove over diluted sgRNAs. Default: 20.

**hit_gene_calling**: method to call hit gene associated with particular phenotype. **Two options: 'position' or 'all', 'position' method is more recommended.** For details, see our paper.

**gene_sgRNA_position**: name of the sgRNA position file (see above, Step 2, File 3).

**Operon_gene_List**: name of the operon file (see above, Step 2, File 4). **This one is optional, if you do not need it, just leave it blank, note that do not delete the tab when leaving this parameter blank**

Below is **an example configure file with default parameters**.

parameter|value
---------|-----
prefix|Ilovemicrobe
fastqpath|example_data
fastq|dCas9R1.fq.gz,dCas9R2.fq.gz,NCR1.fq.gz,NCR2.fq.gz,plasmid.fq.gz
forward_prefixseq|GCAC
forward_suffixseq|GTTT
sample-label|dCas9R1,dCas9R2,NCR1,NCR2,plasmid
sgrna-len|20
list-seq|example_library.csv
experiment_configure|example_experiment_configure.txt
name_configure|example_naming_configure.txt
control_setting|NC
FDR_threshold|0.05
ReadsThreshold|20
hit_gene_calling|position
gene_sgRNA_position|example_coding_region_position.txt
Operon_gene_List|example_operon.txt

After Step 2 and 3, check your working directory. It should looks like below:
[here](./image/files_prepared_before_data_processing.png)

### Step 4：Run the pipeline
Open the command line window (for example, terminal in Macbook), cd to the working directory and run the analysis pipeline.
cd path_to_your_working_directory
python CRISPRscreen_main.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 
python CRISPRscreen_main.py example_configure.txt

Check [here](./image/successful_running.png) for the output during a successful running of the abovementioned test.

The program will create an 'error.log' file under the working directory, open this file to check whether anything wrong happens. Generally, no content suggests successful running. Please post your 'error.log' if you cannot figure out the bugs when using this tool.

For a typical Macbook (for example, 2.6 GHz processor and 8 GB memory), the example test can be finalized within 30 minutes. The rate-limiting step is the mapping of the raw NGS data. For a typical Macbook, we expect a processing speed of 20 million reads per hour. For a genome-scale library with 50 k members (10 sgRNAs per gene assuming 5 k genes encoded by a genome), 100-fold coverage leads to 5 million reads per library, thus roughly 4 NGS libraries processed per hour.

## Output files
The output files will be organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory (prefiex_results). We term this subdirectory 'result directory' thereafter.

[your working directory should be like this after running the test](./image/wkd_after_example_running.png)

You can find many sub directories under the result directory.

[your result directory after running the test](./image/resultdir_after_example_running.png)

Below is the description. For the mathematical processing, see our paper. **All .csv flat files use tab as delimiter unless mentioned**

### NGS raw data profile
-------------------------------------------------------------
#### read count of each sgRNA in each library (prefix_count/)
**prefix.countsummary.txt**: basic statistics of the mapping ratio of each NGS library with a header line using tab as delimiter.

File|Label|Reads|Mapped|Synerror|Unknown|Percentage|Zerocounts|GiniIndex
----|-----|-----|------|--------|-------|----------|----------|---------
example_data/plasmid.fq.gz|plasmid|1000000|838484|90356|71160|0.8385|1734|0.2167
...|...|...|...|...|...|...|...|...

Reads denote the number of reads in the raw data. Mapped denotes number of reads mapping perfectly to one member of the synthetic sgRNA library. Synerror refers to those reads with one indel mutation or more mismatch mutations. Unknown refers to those reads where no forward_prefixseq or forward_suffixseq can be identified. Percentage is the mapping ratio. Zerocount refers to sgRNA number in the *in silico* library without any corresponding read detected. GiniIndex is a metric reflecting the member abundance uniformity in a library. Bigger Gini index indicates more biased distribution with over- represented or diluted members. Generally, more stringent the selective condition is, bigger Gini index we can expect. 

**prefix.count.txt**: raw read count for each sgRNA in the *in silico* library **before normalization**.

sgRNA|dCas9R1|dCas9R2|NCR1|NCR2|plasmid
-----|-------|-------|----|----|-------
gspKb3332_817|12|11|11|9|8
...|...|...|...|...|...

**prefix.normalizeCount.txt**: read count for each sgRNA in the *in silico* library **after normalization of sequencing depth** (for details, see our paper). **This dataset is used for following data processing.**

sgRNA|Gene|dCas9R1|dCas9R2|NCR1|NCR2|plasmid
-----|----|-------|-------|----|----|-------
gspKb3332_817|gspK|12.07|11.08|10.95|8.95|7.98
...|...|...|...|...|...|...

[**prefix_Libray_Gini_Score.png**](./image/all_Libray_Gini_Score.png): schematic of Gini index for each library.

### sgRNA level statistics
-------------------------------------------------------------
#### removed sgRNAs (removed.sgRNA/)
We remove the over diluted sgRNAs with read count less than one threshold ('ReadsThreshold' described in the configure file part).

 1. **prefix.removed.sgRNA.txt**: a simple list flat file with one sgRNA each line

============================================================
#### biological replicate agreement (replicate_consistence/)
Files under this dierectory is a minotoring panel for biological replicate agreement. The replicate information is encoded by the experiment design file (Step 2, File 5). Generally, for N experiments with 2 replicates each, the program produces N scatter plots and N flat files to describe the consistence between replicates for each experiment. One summarizing flat file for all conditions is also given.

 1. **prefix_replicates_reads_statistics.txt**: the summary file
 
 condition|pearson correlation coefficient|P value
 ---------|-------------------------------|-------
 stress1|0.868|0.0
 control1|0.832|0.0
 ...|...|...

 2. **prefix_oneexperiment_replicates.txt**: two-column flat file of read count of each sgRNA in two replicates. For example, NCR1 and NCR2 (in our test data) are two replicates for one experiment.
 
 sgRNA|NCR1_abundance|NCR1_reads|NCR1_abundance_vs_initial|NCR2_abundance|NCR2_reads|NCR2_abundance_vs_initial
 -----|--------------|----------|-------------------------|--------------|----------|-------------------------
 gspKb3332_817|-16.22|10.95|0.46|-16.51|8.95|0.16
 ...|...|...|...|...|...|...

 3. [**prefix_oneexperiment_replicates.png**](./image/all_control1_replicates.png): schematic of replicate agreement of one particular experiment.

============================================================
#### sgRNA read count, abundance change, fitness score, etc (prefix_sgRNA_statistics/)
This directory stores all dataset about sgRNA metrics. It is generally organized at three levels (three sub directories):
 1. **Library level** (information of one NGS library). N files, N = number of rows in experiment design file.
 
 sgRNA|gene|plasmid_Log2_abundnace|plasmid_reads|plasmid_Log2_abundnace_vs_initial
 -----|----|----------------------|-------------|---------------------------------
 gspKb3332_817|gspK|-16.68|7.98|0.0
 ...|...|...|...|...|

 2. **Condition level** (information of one experiement (average of two replicate NGS library)). N files, N = number of columns in experiment design file.
 
 sgRNA|gene|control1_Log2_abundnace|control1_reads|control1_Log2_abundnace_vs_initial|control1_relative_deviation
 -----|----|----------------------|---------------|----------------------------------|---------------------------
 gspKb3332_817|gspK|-16.37|9.90|0.31|0.10
 ...|...|...|...|...|...

 3. **Phenotype level (only need to focus on this level for simplicity)** (information of one phenotype (selective condition normalized by the control condition), under combined_condition_level directory). N files, N = number of 'stress'(selective) conditions in experiment design file.
 
 relative_abundnace_change = Log<sub>2</sub> (read count selective condition / read count control condition)
 
 normalized_change (**it is used as sgRNA fitness score**) = relative_abundnace_change - median relative_abundnace_change of NC sgRNAs
 
 Zscore = normalized_change / sigma of NC sgRNA normalized_change normal distribution
 
 Quality: it is tagged as 'Good' if the averaged read count in control condition is above the threshold ('ReadsThreshold' described in the configure file). **Only 'Good' sgRNAs are used in gene level calculation**.
  
 sgRNA|gene|relative_abundnace_change|normalized_change|Zscore|Quality
 -----|----|-------------------------|-----------------|------|-------
 gspKb3332_817|gspK|0.22|-0.03|-0.04|Good
 ...|...|...|...|...

============================================================
#### NC sgRNA distribution (NCsgRNA_ND/)
Theoretically, fitness socre (log2 abundance change) of NC sgRNA should follow a normal distribution. We hence use a normal distribution to fit NC sgRNA fitness score data.

 1. **prefix_NCsgRNA_ND.txt**: normal distribution of NC sgRNA relative abundance changes (before normalization by median of NC sgRNA relative abundance change, referring to 'relative_abundnace_change' column in phenotype level sgRNA statistics)
 
 condition|median|mean|stdev
 ---------|------|----|-----
 phenotype1|0.25|0.16|0.71
 phenotype2|0.25|0.16|0.71
 ...|...|...|...

 2. **prefix_NCsgRNA_normalized_ND.txt**: normal distribution of NC sgRNA fitness scores (after normalization by median of NC sgRNA relative abundance change, referring to 'normalized_change' column in phenotype level sgRNA statistics)
 
 condition|median|mean|stdev
 ---------|------|----|-----
 phenotype1|0.0|-0.09|0.71
 phenotype2|0.0|-0.09|0.71
 ...|...|...|...

 3. [**prefix_phenotype_NCsgRNAND.png**](./image/all_essential_NCsgRNAND.png): schematic of NC sgRNA fitness score distribution.


### gene level statistics: 
-------------------------------------------------------------
#### FPR-score curve (prefix_quasigeneFPR/)
We use a NC sgRNA derived 'quasi' gene simulation approach (score approach thereafter, we use this method in our paper) (score = |gene fitness| * -Log<sub>10</sub>Pvalue_MWUtest) to determine the false positive rate (*FPR*) for each gene-phenotype association. Hence, for each studied phenotype, the program give 15 simulated *FPR*-score curves with 1 ~ 15 sgRNAs per quasi gene, respectively. Thus, 15 files describing these curves and [one figure file](./image/all_essential_quasigeneFPR.png) are in this sub directory.

============================================================
#### P value-Q value curve (prefix_Pvalue_Qvalue/)
We use a Storey-Tibshirani approach (PNAS 2003) to convert *FPR* into *Q* values. We also use a simple student t test method (sgRNA for one gene vs NC sgRNAs) (t test approach thereafter) to calculate another *P* value and convert it into *Q* values.

 1. **Qvalue_scoreFPR.txt**: score approach derived *P-Q* value curve for all phenotypes.
 
 phenotype1|phenotype2|...|Qvalue
 ----------|----------|---|------
 0.00015|...|...|0.001
 0.00084|...|...|0.005
 ...|...|...|...
 0.03180|...|...|0.1
 
  2. **Qvalue_TtestPvalue.txt**: t test approach derived *P-Q* value curve for all phenotypes, similar to above.
  
  3. Distributions of *FPR* (or t test *P*) for all genes (2N [figures](./image/all_essential_Ttest_pValue.png), N = number of phenotypes).
  
  4. Comparison of *P-Q* curves from two approaches (N [figures](./image/Ilovemicrobe_example_Pvalue_Qvalue.png), N = number of phenotypes).
 
============================================================
#### gene fitness, statistical significance, etc (prefix_gene_statistics/)
This directory stores all dataset about gene metrics. N files are generated, corresponding to N studied phenotypes.

sgRNAnumber: number of sgRNAs used to produce the metrics for this gene (see 'position' approach in our paper)

MedianRAC: gene fintess

MedianZ: gene fitness normalized by the sigma of NC sgRNA fitness normal distribution

-Log10Pvalue_MWUtest: MWU test of sgRNA (of one gene) fitness vs. NC sgRNA fitness

FPRvalue: score approach derived *FPR*; FDRvalue: *Q* value derived from *FPR*;

-Log10Pvalue_Ttest: student t test of sgRNA (of one gene) fitness vs. NC sgRNA fitness

Qvalue_Ttest: *Q* value derived from t test *P* value

gene|sgRNAnumber|MedianRAC|MedianZ|negative Log10Pvalue_MWUtest|FDRvalue|FPRvalue|negative Log10Pvalue_Ttest|Qvalue_Ttest
----|-----------|---------|-------|--------------------|--------|--------|------------------|------------
gspK|11|-0.19|-0.27|0.76|0.30|0.18|0.74|0.36
...|...|...|...|...|...|...|...|...

============================================================
#### FDR-score curve (prefix_quasigeneFDR/)
A [figure](./image/all_essential_quasigeneFDR.png) and a file describing FDR vs. score relation considering sgRNA number per gene profile in the library using a 'quasi' gene simulation approach. It is used to produce the volcano plot.

 phenotype1|phenotype2|...|FDR
 ----------|----------|---|------
 21.91|...|...|0.001
 4.36|...|...|0.005
 ...|...|...|...
 0.72|...|...|0.1

### operon level statistics (reorganize gene data according to the operon file)
-------------------------------------------------------------
#### reorganize gene fitness and statistical significance according to operon structures (prefix_operon_statistics/)
3N files corresponding to N phenotypes. 

**prefix_phenotype_RemoveGene_statistics.txt**: flat file describing removed genes due to the availability of sgRNAs

**prefix_phenotype_RemoveOperon_statistics.txt**: flat file describing removed operons due to the availability of genes

**prefix_phenotype_operon_statistics.txt**: reorganized gene fitness and statistical signifcance dataset. Each line refers to one operon. Different types of metrics are seperated by tab. The same type of metrics corresponding to multiple genes in one operon are seperated by comma. For simplicity, only part of gene metrics are included here.

one phenotype

gene|sgRNAnumber|MedianRAC|MedianZ|negative Log10Pvalue_MWUtest|FDRvalue|FPRvalue
----|-----------|---------|-------|--------------------|--------|--------|
mtlD,mtlR|8,7|-0.110,-0.248|-0.154,-0.347|0.510,0.724|0.472,0.300|0.410,0.184
...|...|...|...|...|...|...|
