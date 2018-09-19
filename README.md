# CRISPRi-functional-genomics-in-prokaryotes

This repository presents the python script collection for CRISPR interference based functional genomics study in prokaryotes, firstly demonstrated in *Escherichia coli* (for the experiment description and result presentation, see https://doi.org/10.1101/129668). This method enables quantitative measurement of fitness score for thousands of genes in multiple conditions, and it works theoretically in any prokaryotic organism. We demonstrated in our paper (Wang T, Guan C, Guo J, Liu B, Wu Y, Xie Z, Zhang C, Xing X-H (2018) Pooled CRISPR interference screening enables genome-scale functional genomics study in bacteria with superior performance. Nat Commun 9:2475 . doi: 10.1038/s41467-018-04899-x) that CRISPR interference exhibits superior performance over other benchmark methods in microbiology research, such as Tn-seq. Hence, we believe it is potent to be a powerful tool for us to better understand the genomics and physiology of microbes. Considering the need of extensive programming and computation to use this method, which is a major hurdle for experimental biologists with limited programming skills, we develop this software package to facilitate the usage of this method by more people.

The python package basically contains two subpackages, one for the sgRNA (N20) design at either customized manner (given interested gene fasta file) or at genome level (given the gene fasta file encoded by particular genome), while another for the NGS data processing and result visualization after screening experiment. For the details about package usage, see the README file within each subpackage depoisited together with the scripts.

[Schematic for the structure of this software package](./package_framework.png)

For any request or question, please contact wtm0217@gmail.com or guanchangge@163.com.

It is also ok to post your concern directly at this GitHub site. This software package is actively maintained and upgraded.

Wish this tool to benefit your work in microbiology study.
