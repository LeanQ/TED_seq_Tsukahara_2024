#!/bin/sh

##Coded by Basile Leduque (basile.leduque@universite-paris-saclay.fr)

######################################################################################
### Input/args ###
######################################################################################

#####input fastq files (paired) are expected to be Gzip compressed and with the following format: name-of-the-sample_1.fq.qz and name-of-the-sample_2.fq.qz. If this is not the case, please modify the first Bowtie2 instance to match the correct file formar.

inname=$1

for inname in filename #change filename
do
echo $inname

inDir=~/HN00212327_TED9 #/path/to/input/fastq/files 
outDir=./"${inname}"_out #/path/to/output/fastq/files ####Modify this to match output folder
refDir=~/TED-seq/ref #/path/to/TAIR10/index
TEseqDir=~/TED-seq/ref/copia93 #/path/to/TAIR10/extremity/index

CORES=4 ##number of CPUs for multithreading. Typically we use 4.
mkdir -p $outDir  

maxcov=1000000  ## here it can be indicated the maximun number of reads requiered to detected an insertion
mincov=50  ## here it can be indicated the minimum number of reads requiered to detected germline insertions. In our experiments, 50 worked well to detect germline insertions. To detect somatic insertions, the pipepeline ignores this parameter
libsize=151  ## here it can be indicated the library size.


#######################################
### 1- Mapping against TE sequence  ###
#######################################

bowtie2-build $TEseqDir/target_TE_sequence_extremity.fa $TEseqDir/target_TE_sequence_extremity

bowtie2 -x $TEseqDir/target_TE_sequence_extremity -1 $inDir/"${inname}"_1.fastq.gz -2 $inDir/"${inname}"_2.fastq.gz -S $inDir/"${inname}".sam --local --very-sensitive

samtools view --threads $CORES -F 1024 -bS $inDir/"${inname}".sam  > $inDir/"${inname}".bam
samtools sort $inDir/"${inname}".bam -o $inDir/"${inname}"-sort.bam
rm $inDir/"${inname}".sam
rm $inDir/"${inname}".bam


########################################
### 2- Extracting informative reads  ###
########################################


### extract pair-end reads in which only one mate mapped ###
samtools view --threads $CORES -F 8 $inDir/"${inname}"-sort.bam  |  awk '{print $1}' | sort -u > $outDir/"${inname}"_reads.name.disc
samtools view -f 64 -u $inDir/"${inname}"-sort.bam > $inDir/"${inname}"-disc-first.bam 
samtools view -f 128 -u $inDir/"${inname}"-sort.bam > $inDir/"${inname}"-disc-second.bam 

singularity exec -B /usr/lib/locale/locale-archive:/usr/lib/locale/locale-archive -B /usr/bin/locale:/usr/bin/locale /usr/local/biotools/p/picard:2.27.5--hdfd78af_0 picard -Xmx3g FilterSamReads INPUT=$inDir/"${inname}"-disc-first.bam FILTER=includeReadList READ_LIST_FILE=$outDir/"${inname}"_reads.name.disc OUTPUT=$outDir/"${inname}"selected-disc-first.sam 
singularity exec -B /usr/lib/locale/locale-archive:/usr/lib/locale/locale-archive -B /usr/bin/locale:/usr/bin/locale /usr/local/biotools/p/picard:2.27.5--hdfd78af_0 picard -Xmx3g  FilterSamReads INPUT=$inDir/"${inname}"-disc-second.bam FILTER=includeReadList READ_LIST_FILE=$outDir/"${inname}"_reads.name.disc OUTPUT=$outDir/"${inname}"selected-disc-second.sam 



cat $outDir/"${inname}"selected-disc-first.sam  | awk '$1!~/^@/ {print $1"\t"$10"\t"$11}' | sort -u -k1,1 -k2,2 | awk '{print "@"$1"|1\n"$2"\n+\n"$3}' > $outDir/"${inname}"-disc.fastq 
cat $outDir/"${inname}"selected-disc-second.sam | awk '$1!~/^@/ {print $1"\t"$10"\t"$11}' | sort -u -k1,1 -k2,2 | awk '{print "@"$1"|2\n"$2"\n+\n"$3}' >> $outDir/"${inname}"-disc.fastq 
    
    
### extract reads that map discordantly ###

samtools view --threads $CORES -hF 4 $inDir/"${inname}"-sort.bam | awk -v l=$LS '(($4-$8)>(l*10) || ($8-$4)>(l*10) || $7!="=") || $1~/@HD/ || $1~/@SQ/ || $1~/@PG/' >  $outDir/"${inname}"_clip.sam 
samtools sort --threads $CORES -n  $outDir/"${inname}"_clip.sam | samtools fastq --threads $CORES -F 256 /dev/stdin -1  $outDir/"${inname}"_clip_reads_1.fq -2  $outDir/"${inname}"_clip_reads_2.fq
cat   $outDir/"${inname}"_clip_reads_1.fq  $outDir/"${inname}"_clip_reads_2.fq  $outDir/"${inname}"-disc.fastq > $outDir/"${inname}"_clip_disc_reads.fq
rm $outDir/"${inname}"_clip_reads_1.fq 
rm $outDir/"${inname}"_clip_reads_2.fq 
rm  $outDir/"${inname}"-disc.fastq 


#########################################################
### 3- mapping informative reads on reference genome  ###
#########################################################


bowtie2 -x $refDir/Col-CEN_v1.2 -U  $outDir/"${inname}"_clip_disc_reads.fq -S $outDir/"${inname}"_clip_disc-local.sam --local --very-sensitive -p $CORES
samtools view --threads $CORES -bS $outDir/"${inname}"_clip_disc-local.sam > $outDir/"${inname}"_clip_disc-local.bam
samtools sort --threads $CORES $outDir/"${inname}"_clip_disc-local.bam > $outDir/"${inname}"_clip_disc-local.sorted.bam
samtools index $outDir/"${inname}"_clip_disc-local.sorted.bam



##############################################################
### 4- Detecting somatic insertions 
##############################################################
## from the  "${inname}"_clip_disc-local.sorted.bam generate a bed file ($outDir/"${inname}"_somatic_ins.bed) containing the localisation of every insertion (Somatic + germline)


Qscore=5  ## Quality score threshold,  5 call insertion supported by uniquely mapped reads 
samtools view -H $outDir/"${inname}"_clip_disc-local.sorted.bam > $outDir/"${inname}"_somatic.sorted.sam
samtools view  --threads $CORES -q $Qscore  $outDir/"${inname}"_clip_disc-local.sorted.bam | awk  '$6~/S/ && $1~"\\|1$" { print $0}' >> $outDir/"${inname}"_somatic.sorted.sam #&& $1~"\\|1$"を追加
samtools view  --threads $CORES -bS  $outDir/"${inname}"_somatic.sorted.sam > $outDir/"${inname}"_somatic.sorted.bam

bamToBed -i $outDir/"${inname}"_somatic.sorted.bam -split > $outDir/"${inname}"_somatic_bamtobed.bed

awk '{ print $1"\t"$2"\t"$3"\t"$6}' $outDir/"${inname}"_somatic_bamtobed.bed | awk '!seen[$0]++' > $outDir/"${inname}"_uniq_somatic_bamtobed.bed

awk  '$4=="+" { print $1"\t"$3-1"\t"$3} '  $outDir/"${inname}"_uniq_somatic_bamtobed.bed | awk '!seen[$0]++' > $outDir/"${inname}"_somatic_ins.bed
awk  '$4=="-" { print $1"\t"$2"\t"$2+1} '  $outDir/"${inname}"_uniq_somatic_bamtobed.bed | awk '!seen[$0]++' >> $outDir/"${inname}"_somatic_ins.bed
sort  -u -k1,1 -k2,2n $outDir/"${inname}"_somatic_ins.bed  > $outDir/"${inname}"_somatic_ins_sorted.bed

rm $outDir/"${inname}"_somatic.sorted.bam
rm $outDir/"${inname}"_uniq_somatic_bamtobed.bed
rm $outDir/"${inname}"_somatic_bamtobed.bed
rm $outDir/"${inname}"_somatic.sorted.sam

##############
bedtools intersect  -v  -a $outDir/"${inname}"_somatic_ins_sorted.bed  -b $TEseqDir/targeted_TE_sequences.bed > $outDir/"${inname}"_somatic_ins_rmevade_sorted.bed
########################
bedtools coverage -a $refDir/Col-CEN_10kb_windows.bed -b $outDir/"${inname}"_somatic_ins_rmevade_sorted.bed > $outDir/"${inname}"_somatic_ins_10kb_windows.bed

##############################################################
### 5- generating genome-wide profiles of insertions in 10Kbp windows 
##############################################################

#bedtools makewindows -g Col-CEN_v1.2_file.txt -w 90000 -s 10000 > $refDir/Col-CEN_10kb_sliding_windows_size9.bed
bedtools coverage -a $refDir/Col-CEN_10kb_sliding_windows_size9.bed -b $outDir/"${inname}"_somatic_ins_rmevade_sorted.bed > $outDir/"${inname}"_somatic_ins_rmevade_size9.bed
awk  '{ print $1"\t"$2"\t"$3"\t"$4/9"\t"$5"\t"$6"\t"$7} ' $outDir/"${inname}"_somatic_ins_rmevade_size9.bed > $outDir/"${inname}"_somatic_ins_rmevade_10kb_sliding_size9.bed

rm $outDir/"${inname}"_somatic_ins_rmevade_size9.bed
rm $outDir/"${inname}"_somatic.sorted.bam
rm $outDir/"${inname}"_somatic.sorted.sam


done
#
