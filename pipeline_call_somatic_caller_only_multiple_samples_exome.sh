#!/bin/bash

#pipeline to analyse paired samples for somatic and mosacic variants
#Alan Pittman September 2018

#requires input blood and tumour_options ~ please specify

#usage as follows:
#sh Pipeline_mozaic -n normal/blood-sample -t tumor/biopsy-sample

########################################################################################
########################################################################################
#resources

BWAindex=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/human_g1k_v37.fasta
java=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/java/jre1.8.0_171/bin/java
picard=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/picard-2.815/picard.jar
gatk=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar
ExomeTarget=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/BroadExACExomeIntervlas.bed
NRAS=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/NRAS.bed
refknownsitesSNPS=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/common_all.vcf

########################################################################################
#limiting analysis to particular genomic intervals with -L

project=$1
InBamsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/Exomes/Aligned/$project #targeted resequencing  
mutect2Out=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/raw_somatic_caller_only


sleep 1
mkdir $mutect2Out/$project
samples=`ls $InBamsDir`

echo $samples

sleep 5

for sample in $samples; do

inBAM=$InBamsDir/$sample/${sample}_sorted_unique_recalibrated.bam

echo $inBAM

#check your exome definition here

$java -jar $gatk Mutect2 \
	-R $BWAindex \
	-I $inBAM \
	-tumor $sample \
	-L $ExomeTarget \
	-O $mutect2Out/$project/${sample}_raw_muTect2.vcf \
	-bamout $mutect2Out/$project/${sample}_.bam 
	
done	
	
exit
exit



