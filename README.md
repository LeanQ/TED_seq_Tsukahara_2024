# TED-seq Pipeline to detect somatic insertions used in Tsukahara et al 2024

This pipeline was implemented to discover somatic TE insertions based on TED-seq sequencing data.


### Dependencies

SPLITREADER requires the following softwares:

* SAMtools (v1.2 or higher) (http://samtools.sourceforge.net/)
* Bowtie2 (v2.2.9 or higher) https://sourceforge.net/projects/bowtie-bio/files/bowtie2/)
* Picard tools (> Java 1.8) (https://broadinstitute.github.io/picard/)
* bedtools (v2.20.1 or higher) (https://bedtools.readthedocs.io/en/latest/content/installation.html)

## Main steps

The pipeline consists in four steps: 

1- Forced mapping to targeted TE extremeties captured by TED-seq "target_TE_sequence_extremity"

2- Extracting informative reads over targeted TE extremeties

3- Mapping informative reads on reference genome 

4- Detecting somatic insertions

5- Generation of genome-wide profiles of somatic insertions in 10Kbp windows 


## Input data

This pipeline requieres as an input:

1. Pair-end fastq files from a TED-seq experiment
   
2. Targeted TE extremity sequence and bowtie2 indexes
   
3. Reference genome sequence and and bowtie2 indexes
  
4. A bed file containing 10Kbp windows covering the whole referene genome. 

