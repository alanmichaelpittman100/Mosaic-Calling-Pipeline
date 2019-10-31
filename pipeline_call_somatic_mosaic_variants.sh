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

inBAMsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_BAMs/
PoNVCF=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_vcfs/
PoNjointVCFDIR=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_joint_vcfs/

mutect2Out=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/raw_somatic_calls
InBamsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/Exomes/Aligned/*


PoNjointVCF="PoN.vcf.gz" #specify what panel of normal vcf files you want to use

while getopts n:t: option
do
        case "${option}"
        in

				n) myIDnormal=${OPTARG};;
				t) myIDtumor=${OPTARG};;

  esac

done

#now lets do the paired somatic calling:

normalBAM=$InBamsDir/$myIDnormal/${myIDnormal}_sorted_unique_recalibrated.bam
tumorBAM=$InBamsDir/$myIDtumor/${myIDtumor}_sorted_unique_recalibrated.bam


$java -jar $gatk Mutect2 \
	-R $BWAindex \
	-I $normalBAM \
	-I $tumorBAM \
	-normal $myIDnormal \
	-tumor $myIDtumor \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--af-of-alleles-not-in-resource 0.0000025 \
	-pon ${PoNjointVCFDIR}/$PoNjointVCF \
	-L $ExomeTarget \
	-O $mutect2Out/${myIDnormal}_${myIDtumor}_raw_muTect2.vcf \
	-bamout $mutect2Out/${myIDnormal}_${myIDtumor}_.bam

	
exit
exit

