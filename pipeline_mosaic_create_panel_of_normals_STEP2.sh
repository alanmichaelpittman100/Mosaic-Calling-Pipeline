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
PoNjointVCF=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_joint_vcfs

PoNsamples="M0008
M0009
M0010
M0012
M0014
M0017
M0018
M0023
M0027
S2459
S2469
S2473"


rm ${PoNjointVCF}/gVCFs.txt

for sample in $PoNsamples; do

	echo "--vcfs $PoNVCF/${sample}/${sample}.vcf.gz" >> ${PoNjointVCF}/gVCFs.txt
	
done

#now lets combine all the vcfs into 1

$java -jar $gatk \
 CreateSomaticPanelOfNormals \
 `cat ${PoNjointVCF}/gVCFs.txt` \
 -O ${PoNjointVCF}/PoN.vcf.gz
	
exit
exit
exit
