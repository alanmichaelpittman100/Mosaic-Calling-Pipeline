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
refknownsitesSNPS=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/common_all.vcf

########################################################################################
########################################################################################

project=$1

InBamsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/Exomes/Aligned/$project
mutect2Out=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/raw_somatic_caller_only/$project
filtermutect2out=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/filtered_somatic_calls/$project

sleep 1
mkdir $filtermutect2out
samples=`ls $InBamsDir`

echo $samples

sleep 5

for sample in $samples; do

$java -jar $gatk FilterMutectCalls \
	-V $mutect2Out/${sample}_raw_muTect2.vcf \
	-O $filtermutect2out/${sample}_filtered_muTect2.vcf
done

exit
