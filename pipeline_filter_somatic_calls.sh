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

PoNjointVCF="PoN.vcf.gz" #specify what panel of normal vcf files you want to use

########################################################################################

inBAMsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_BAMs/
PoNVCF=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_vcfs/
PoNjointVCFDIR=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/panel_of_normals_joint_vcfs/

mutect2Out=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/raw_somatic_calls
InBamsDir=/homes/athosnew/Genetics_Centre_Bioinformatics/Exomes/Aligned/*


mutect2Filtered=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/filtered_somatic_calls

contaminationtablesDir=/homes/athosnew/Genetics_Centre_Bioinformatics/mosaic_calling/contamination_tables

commonVariants=/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/hapmap_3.3.b37.vcf



while getopts n:t: option
do
        case "${option}"
        in

				n) myIDnormal=${OPTARG};;
				t) myIDtumor=${OPTARG};;

  esac

done

tumorBAM=$InBamsDir/$myIDtumor/${myIDtumor}_sorted_unique_recalibrated.bam

#STEP1 get pileup summaries:

$java -jar $gatk GetPileupSummaries \
	-I $tumorBAM \
	-V $commonVariants \
	-O $contaminationtablesDir/${myIDtumor}_getpileupsummaries.table

#STEP2 calculate contamination:

$java -jar $gatk CalculateContamination \
	-I $contaminationtablesDir/${myIDtumor}_getpileupsummaries.table \
	-O $contaminationtablesDir/${myIDtumor}_calculatecontamination.table

#STEP3 Filter MuTect2 calls:

$java -jar $gatk FilterMutectCalls \
	-V $mutect2Out/${myIDnormal}_${myIDtumor}_raw_muTect2.vcf \
	--contamination-table $contaminationtablesDir/${myIDtumor}_calculatecontamination.table \
	-O $mutect2Filtered/${myIDnormal}_${myIDtumor}_filtered_muTect2.vcf

exit
exit


