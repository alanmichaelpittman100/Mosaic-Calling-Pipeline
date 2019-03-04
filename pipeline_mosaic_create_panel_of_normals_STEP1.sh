#!/bin/bash

#pipeline to create a panel of normals for muTect2
#Alan Pittman October 2018

########################################################################################
#resources

BWAindex=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/human_g1k_v37.fasta
java=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/java/jre1.8.0_171/bin/java
picard=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/picard-2.815/picard.jar
gatk=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar

ExomeTarget=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/BroadExACExomeIntervlas.bed

########################################################################################

inBAMsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_BAMs/
PoNVCF=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_vcfs/
PoNjointVCF=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_joint_vcfs/

PoNsamples="MyNormalSample1 MyNormalSample2 MyNormalSample3 MyNormalSample4" # works best with 100 control exomes

for sample in $PoNsamples; do

	mkdir $PoNVCF/${sample}/

	$java -jar $gatk Mutect2 \
	-R $BWAindex \
	-I ${inBAMsDir}/${sample}/${sample}_sorted_unique_recalibrated.bam \
	-tumor ${sample} \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	-L $ExomeTarget \
	-O $PoNVCF/${sample}/${sample}.vcf.gz

done


exit
exit
