# Mosaic-Calling-Pipeline
Protocol for identification of somatic/mosaic variants in blood-biopsy pairs

GATK MuTect2 workflow

##Requirements

human_g1k_v37.fasta

java

picardtools (picard-2.815 tested here) 

gatk (gatk-4.0.4.0 tested here)

gatk resource bundle (https://github.com/bahlolab/bioinfotools/blob/master/GATK/resource_bundle.md)

common_all.vcf (ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/)


##step 1: create panel of normals

pipeline_mosaic_create_panel_of_normals_STEP1.sh

pipeline_mosaic_create_panel_of_normals_STEP2.sh

##step 2: call mosaic variants

pipeline_call_somatic_mosaic_variants.sh

##step 3: filter somatic calls 

pipeline_filter_somatic_calls.sh

