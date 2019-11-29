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


#test samples S1575 and S1574


$java -jar $gatk \
 CreateSomaticPanelOfNormals \
 `cat ${PoNjointVCF}/gVCFs.txt` \
 -O ${PoNjointVCF}/PoN.vcf.gz




gatk --java-options "-Xmx2g" Mutect2 \
-R hg38/Homo_sapiens_assembly38.fasta \
-I tumor.bam \
-I normal.bam \
-tumor HCC1143_tumor \
-normal HCC1143_normal \
-pon resources/chr17_pon.vcf.gz \
--germline-resource resources/chr17_af-only-gnomad_grch38.vcf.gz \
--af-of-alleles-not-in-resource 0.0000025 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 1_somatic_m2.vcf.gz \
-bamout 2_tumor_normal_m2.bam 











##useful stuff;


########################################################################################
#resources:

novoalign="/data/kronos/NGS_Software/novocraft_v3/novoalign"
novoalignArguments="-c 8 --rOQ --hdrhd 3 -H -k -o Soft -t 320 -F ILM1.8"
indexedgenome="/data/kronos/NGS_Reference/novoalign/human_g1k_v37.fasta.k15.s2.novoindex"
samtools="/data/kronos/NGS_Software/samtools-0.1.18/samtools"
picard="/data/kronos/NGS_Software/picard-tools-1.75"
gatk="/data/kronos/NGS_Software/GATK_v3_3/GenomeAnalysisTK.jar"
gatk5="/data/kronos/NGS_Software/GATK_v3_5/GenomeAnalysisTK.jar"
genomeFASTA="/data/kronos/NGS_Reference/fasta/human_g1k_v37.fasta"
genomeFASTAI="/data/kronos/NGS_Reference/fasta/human_g1k_v37.fastai"
JAVA="/data/kronos/General_Software/jre1.7.0_67/bin/java"
knownINDELS="/data/kronos/NGS_Reference/GATK_refFiles/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf"
samtoolsMEM="5000000000" 
picardArguments="TMP_DIR=/data/kronos/temp ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE"
INTERVALS_AgilentV6="/data/kronos/NGS_Reference/GATK_refFiles/SureSelect_V6_UTR.bed"
VARSCANOPT="--min-avg-qual 20 --tumor-purity 0.30 --normal-purity 0.99 --min-coverage-tumor 12 --strand-filter 1 --min-var-freq 0.15 --min-strands2 2"  #tinker around with these :-) 

########################################################################################

oFolder="/array/apittman/VK_mosaic_Analysis/Results_output/AVM"

while getopts n:t: option
do
        case "${option}"
        in

				n) myIDblood=${OPTARG};;
				t) myIDtumor=${OPTARG};;

  esac

done

mkdir $oFolder/${myIDblood}_${myIDtumor}; echo "making output directories"

mkdir $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS

mkdir $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants
mkdir $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants

mkdir $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2
mkdir $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2


sleep 2

iFolder="/array/apittman/VK_mosaic_Analysis/BAM_repository"

echo "#!/bin/bash" > $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP1 - Joint realignment around INDELS
#########################################################################################

echo "java -jar $gatk5 -T RealignerTargetCreator -nt 8 -R $genomeFASTA -L $INTERVALS_AgilentV6 --known $knownINDELS -I ${iFolder}/${myIDblood}/${myIDblood}_sorted_unique_realigned.bam -I ${iFolder}/${myIDtumor}/${myIDtumor}_sorted_unique_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_sorted_unique_realigned.bam.list" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "java -jar $gatk5 -T IndelRealigner -targetIntervals $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_sorted_unique_realigned.bam.list -R $genomeFASTA --knownAlleles $knownINDELS -I ${iFolder}/${myIDblood}/${myIDblood}_sorted_unique_realigned.bam -I ${iFolder}/${myIDtumor}/${myIDtumor}_sorted_unique_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_joint_realigned.bam" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#Dispence individual bam files:

#GATK print reads Command:

echo "java -jar $gatk5 -T PrintReads -R $genomeFASTA -nct 8 --sample_name ${myIDblood} -I $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_joint_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_realigned.bam" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "java -jar $gatk5 -T PrintReads -R $genomeFASTA -nct 8 --sample_name ${myIDtumor} -I $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_joint_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDtumor}_realigned.bam" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

echo "rm $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_joint_realigned.bam" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "rm $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_joint_realigned.bai" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "rm $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_${myIDtumor}_sorted_unique_realigned.bam.list" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP2 - call germline variants - HaplotypeCaller
#########################################################################################

echo "java -jar $gatk5 -T HaplotypeCaller -nct 8 -R $genomeFASTA -L $INTERVALS_AgilentV6 -pcrModel CONSERVATIVE -mbq 10 -stand_call_conf 30 -stand_emit_conf 10 -out_mode EMIT_VARIANTS_ONLY -I $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "java -jar /data/kronos/NGS_Software/GATK_v3_5/GenomeAnalysisTK.jar -T VariantFiltration -R $genomeFASTA --variant $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 

echo "/data/kronos/NGS_Software/annovar/convert2annovar.pl -format vcf4old $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf -includeinfo > $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "/data/kronos/NGS_Software/annovar/table_annovar.pl $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput /data/kronos/NGS_Software/annovar_Nov2014/humandb_hg19/ -buildver hg19 -protocol refGene,genomicSuperDups,1000g2012apr_all,esp6500si_all,esp6500si_aa,esp6500_ea,snp129,snp137,cg69,exac02,ljb_all,clinvar_20140211,caddgt20,DMN_Exomes_Freq -remove -otherinfo -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring n_a" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

echo "java -jar $gatk5 -T HaplotypeCaller -nct 8 -R $genomeFASTA -L $INTERVALS_AgilentV6 -pcrModel CONSERVATIVE -mbq 10 -stand_call_conf 30 -stand_emit_conf 10 -out_mode EMIT_VARIANTS_ONLY -I $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDtumor}_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "java -jar $gatk5 -T VariantFiltration -R $genomeFASTA --variant $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf -o $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 

echo "/data/kronos/NGS_Software/annovar/convert2annovar.pl -format vcf4old $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf -includeinfo > $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "/data/kronos/NGS_Software/annovar/table_annovar.pl $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput /data/kronos/NGS_Software/annovar_Nov2014/humandb_hg19/ -buildver hg19 -protocol refGene,genomicSuperDups,1000g2012apr_all,esp6500si_all,esp6500si_aa,esp6500_ea,snp129,snp137,cg69,exac02,ljb_all,clinvar_20140211,caddgt20,DMN_Exomes_Freq -remove -otherinfo -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring n_a" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#Variant Evaluation; compare the concordnace of blood calls to skin calles ; do the samples match ??????
echo "vcftools --vcf $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf --diff $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf --out $oFolder/${myIDblood}_${myIDtumor}/compared_germline_variants_numbers --diff-site --bed /data/kronos/NGS_Reference/GATK_refFiles/INTERVALS_focussed.bed --remove-filtered-all" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP2 - muTECT2 analysis
#########################################################################################

echo "java -jar $gatk5 -T MuTect2 -L $INTERVALS_AgilentV6 --min_base_quality_score 20 -R $genomeFASTA  -I:tumor $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDtumor}_realigned.bam -I:normal $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_realigned.bam -o $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2/${myIDblood}_${myIDtumor}_muTECT2.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
#nct 8 taken out
#notes:
#--normal_panel
#--max_alt_alleles_in_normal_count

echo "/data/kronos/NGS_Software/annovar/convert2annovar.pl -format vcf4old $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2/${myIDblood}_${myIDtumor}_muTECT2.vcf -includeinfo > $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2/${myIDblood}_${myIDtumor}_muTECT2.vcf.avinput" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "/data/kronos/NGS_Software/annovar/table_annovar.pl $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2/${myIDblood}_${myIDtumor}_muTECT2.vcf.avinput /data/kronos/NGS_Software/annovar_Nov2014/humandb_hg19/ -buildver hg19 -protocol refGene,genomicSuperDups,1000g2012apr_all,esp6500si_all,esp6500si_aa,esp6500_ea,snp129,snp137,cg69,exac02,ljb_all,clinvar_20140211,caddgt20,DMN_Exomes_Freq -remove -otherinfo -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring n_a" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP4 - Varsacn 2 analysis
#########################################################################################
echo "/data/kronos/NGS_Software/samtools-1.3.1/samtools mpileup -q 20 -Q 20 -l $INTERVALS_AgilentV6 -f $genomeFASTA $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDblood}_realigned.bam  > $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}.pileup" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "/data/kronos/NGS_Software/samtools-1.3.1/samtools mpileup -q 20 -Q 20 -l $INTERVALS_AgilentV6 -f $genomeFASTA $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Joint_realigned_BAMS/${myIDtumor}_realigned.bam  > $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDtumor}.pileup" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#where does the VCF end up ????????????? in the current directory where you run the script from ???
echo "java -jar /data/kronos/NGS_Software/Varscan2/VarScan.v2.4.1.jar somatic $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}.pileup $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDtumor}.pileup $VARSCANOPT --output-vcf 1" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 
echo "cp $oFolder/${myIDblood}_${myIDtumor}/.vcf $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 
echo "rm -rf $oFolder/${myIDblood}_${myIDtumor}/.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 

#rm $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}.pileup
#rm $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDtumor}.pileup

#select out 'SOMATIC' variants

echo "grep '#' $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}.vcf >> $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_somatic.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 
echo "grep 'SOMATIC' $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}.vcf >> $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_VARSCAN_somatic.vcf" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh 

#annotate VARSCAN2 .vcf

echo "/data/kronos/NGS_Software/annovar/convert2annovar.pl -format vcf4old $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_VARSCAN_somatic.vcf -includeinfo > $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_VARSCAN_somatic.vcf.avinput" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "/data/kronos/NGS_Software/annovar/table_annovar.pl $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_VARSCAN_somatic.vcf.avinput /data/kronos/NGS_Software/annovar_Nov2014/humandb_hg19/ -buildver hg19 -protocol refGene,genomicSuperDups,1000g2012apr_all,esp6500si_all,esp6500si_aa,esp6500_ea,snp129,snp137,cg69,exac02,ljb_all,clinvar_20140211,caddgt20,DMN_Exomes_Freq -remove -otherinfo -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring n_a" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP5 - LoFreq analysis
#########################################################################################

#?????????????????????????????????????????????????????????????????????????????????????????

#STEP6 - Final Processing - variant filtering- make a nice excel file ! -python?
#########################################################################################

#call python script
echo "python /array/apittman/VK_mosaic_Analysis/Final_variant_processing_v2.py $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2/${myIDblood}_${myIDtumor}_muTECT2.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_VARSCAN_somatic.vcf.avinput.hg19_multianno.txt" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP6 - Final Processing - variant filtering- make a nice excel file ! -python?
#########################################################################################


#make final folder of short name files:
echo "mkdir $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

echo "cp $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_germlineVariants/${myIDblood}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/${myIDblood}_GER.txt" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "cp $oFolder/${myIDblood}_${myIDtumor}/${myIDtumor}_germlineVariants/${myIDtumor}.raw.snps.indels.vcf_VARIANT_FILTRATION.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/${myIDtumor}_GER.txt" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "cp $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_muTECT2/${myIDblood}_${myIDtumor}_muTECT2.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/muTECT.txt" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
echo "cp $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_VARSCAN2/${myIDblood}_${myIDtumor}_VARSCAN_somatic.vcf.avinput.hg19_multianno.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/Varscan2.txt" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#call python script
echo "python /array/apittman/VK_mosaic_Analysis/Final_variant_processing_v2.py $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/${myIDblood}_GER.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/${myIDtumor}_GER.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/muTECT.txt $oFolder/${myIDblood}_${myIDtumor}/${myIDblood}_${myIDtumor}_Final_Analysis_files/Varscan2.txt" >> $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh

#STEP - FINAL - submit jobs to SGE on KRONOS
#########################################################################################

sleep 1
echo "
	=======================================
	       Submitting jobs to Kronons
	=======================================
	"
sleep 1

qsub -pe make 6 -cwd $oFolder/${myIDblood}_${myIDtumor}/J_script_pipeline_mozaic_${myIDblood}_${myIDtumor}.sh
				
exit
exit
exit
exit


