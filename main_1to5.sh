#!/bin/bash

## Following instructions to install ensembleCNV, PennCNV.
## Assume the installation folder
PIPELINE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2"
PENNCNV="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/PennCNV-1.0.5/"
FINALREPORT="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_final_report_16022022_Ovine50K_4.txt"
SNPMAP="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_snp_map_16022022_Ovine50K_4.txt"
PHENOTYPE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/Pag_sheep_phenotypes_OneRecordPerParity.txt"
#====================================================================================================
ENSEMBLECNV=${PIPELINE}/04_Run_EnsembleCNV
## Working directory for a ensemble new project
WKDIR=${PIPELINE}/04a_EnsembleCNV_Working_Directory
#====================================================================================================
cd ${PIPELINE}

module load python/3.11.4
module load R/4.3.0

## Define the directory for R libraries
R_LIBS_DIR="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/R_libs/"

## Check if the directory exists, if not, create it
if [ ! -d "$R_LIBS_DIR" ]; then
    mkdir -p "$R_LIBS_DIR"
fi

## Set R_LIBS environment variable
export R_LIBS=$R_LIBS_DIR

## Install necessary R packages if not already installed
Rscript -e "
packages <- c('cowplot', 'data.table', 'dplyr', 'ggplot2', 'gridExtra', 'mclust', 'mixtools', 'modeest', 'optparse', 'pheatmap', 'plyr', 'RColorBrewer', 'Rcpp', 'RcppArmadillo', 'tibble')
install_if_missing <- function(p) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, lib = '$R_LIBS_DIR', repos = 'http://cran.us.r-project.org')
  }
}
invisible(sapply(packages, install_if_missing))
"
#****************************************PART1 Update Input Data****************************************
mkdir -p ${PIPELINE}/01_Data/results/

# update snp
python ${PIPELINE}/01_Data/scripts/step1.update_snp_map_to_3.1.py \
${PIPELINE}/01_Data/data/SNPchimp_result_2721813773.tsv \
${SNPMAP} \
${PIPELINE}/01_Data/results/new_snp_map_16022022_Ovine50K_4.txt

#====================================================================================================
# update final report
python ${PIPELINE}/01_Data/scripts/step2.update_final_report_to_3.1.py \
${PIPELINE}/01_Data/results/new_snp_map_16022022_Ovine50K_4.txt \
${FINALREPORT} \
${PIPELINE}/01_Data/results/new_final_report_16022022_Ovine50K_4.txt

#****************************************PART2 Quality Control****************************************
mkdir -p ${PIPELINE}/02_Quality_Control_by_Plink/results

# generate files for plink
python ${PIPELINE}/02_Quality_Control_by_Plink/scripts/step1.generate_ped_map.py \
  --genotype_file "${PIPELINE}/01_Data/results/new_final_report_16022022_Ovine50K_4.txt" \
  --snp_map_file "${PIPELINE}/01_Data/results/new_snp_map_16022022_Ovine50K_4.txt" \
  --ped_file "${PIPELINE}/02_Quality_Control_by_Plink/results/input_data.ped" \
  --map_file "${PIPELINE}/02_Quality_Control_by_Plink/results/input_data.map"

#====================================================================================================
# run plink
cd ${PIPELINE}/02_Quality_Control_by_Plink/results

bash ${PIPELINE}/02_Quality_Control_by_Plink/scripts/step2.run_plink.sh \
${PIPELINE}/02_Quality_Control_by_Plink/results

#****************************************PART3 Generate Input File for EnsembleCNV****************************************
## Create a new project
cd $ENSEMBLECNV
chmod +x create_new_project.sh

./create_new_project.sh $WKDIR

mkdir -p ${PIPELINE}/03_Generate_Input_Files/results
## copy centromere to working directory
cp ${PIPELINE}/01_Data/data/centromere.txt ${WKDIR}/data/centromere.txt

#====================================================================================================
# generate final report files for ensembleCNV
python ${PIPELINE}/03_Generate_Input_Files/scripts/step2.generate_final_report.py \
--report_path ${PIPELINE}/01_Data/results/new_final_report_16022022_Ovine50K_4.txt \
--map_path ${PIPELINE}/01_Data/results/new_snp_map_16022022_Ovine50K_4.txt \
--output_path ${PIPELINE}/03_Generate_Input_Files/results/no_qc_final_report.txt \
--missing_data_output_path ${PIPELINE}/03_Generate_Input_Files/results/missing_data_report.txt

#====================================================================================================
# clean final report files for ensembleCNV based on plink results
python ${PIPELINE}/03_Generate_Input_Files/scripts/step3.clean_final_report.py \
--final_report_path ${PIPELINE}/03_Generate_Input_Files/results/no_qc_final_report.txt \
--qc_individuals_path ${PIPELINE}/02_Quality_Control_by_Plink/results/plink_results/filtered_data_final_individuals.txt \
--qc_snps_path ${PIPELINE}/02_Quality_Control_by_Plink/results/plink_results/filtered_data_final_snps.txt \
--output_path ${WKDIR}/data/final_report.txt

#====================================================================================================
# generate samples table file for ensembleCNV
Rscript ${PIPELINE}/03_Generate_Input_Files/scripts/step4.generate_samples_table.R \
${WKDIR}/data/final_report.txt \
${WKDIR}/data/Samples_Table.txt

#****************************************PART4 Run EnsembleCNV****************************************

## 1 Initial call =============================================================

### Prepare chromosome-wise LRR and BAF matrices for CNV genotyping -----------

#### (1) Create LRR and BAF (tab delimited) matrices from final report
perl ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/finalreport_matrix_LRR_BAF.pl \
${WKDIR}/data/final_report.txt \
${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF

#### (2) Tansform tab-delimited text file to .rds format for quick loading in R
Rscript ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/transform_from_tab_to_rds.R \
--input ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF \
--output ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \

### Prepare data for individual CNV callers -----------------------------------
#### PennCNV
perl ${WKDIR}/01_initial_call/prepare_IPQ_input_file/finalreport_to_PennCNV.pl \
-prefix ${WKDIR}/01_initial_call/run_PennCNV/data/ \
-suffix .txt \
--tolerate \
${WKDIR}/data/final_report.txt

## run_PennCNV ----------------------------------------------------------------

#### (1) Prepare SNP.pfb
#### compile pfb (population frequency of B allele) file

if [ ! -f "${WKDIR}/01_initial_call/run_PennCNV/data_aux/list_pfb.txt" ]; then
    touch "${WKDIR}/01_initial_call/run_PennCNV/data_aux/list_pfb.txt"
fi
find "${WKDIR}/01_initial_call/run_PennCNV/data" -type f -name "*.txt" > "${WKDIR}/01_initial_call/run_PennCNV/data_aux/list_pfb.txt"

perl ${PENNCNV}/compile_pfb.pl \
-snpposfile ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt \
-listfile ${WKDIR}/01_initial_call/run_PennCNV/data_aux/list_pfb.txt \
-output ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb

#### (2) Run PennCNV for each sample in parallel (through job submitting system on cluster)
Rscript ${WKDIR}/01_initial_call/run_PennCNV/step.2.run.PennCNV.jobs.R \
--penncnv ${PENNCNV} \
--data ${WKDIR}/01_initial_call/run_PennCNV/data \
--wkdir ${WKDIR}/01_initial_call/run_PennCNV/results \
--pfb ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb \
--hmm ${PENNCNV}/lib/hhall.hmm

#### (3) Check job status and resubmit failed jobs
Rscript ${WKDIR}/01_initial_call/run_PennCNV/step.3.check.PennCNV.jobs.R \
--penncnv ${PENNCNV} \
--data ${WKDIR}/01_initial_call/run_PennCNV/data/ \
--wkdir ${WKDIR}/01_initial_call/run_PennCNV/results \
--pfb ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb \
--hmm ${PENNCNV}/lib/hhall.hmm

#### (4) Combine PennCNV results (.rawcnv and .log files) from each sample
perl ${WKDIR}/01_initial_call/run_PennCNV/step.4.combine.PennCNV.res.pl \
--in_dir ${WKDIR}/01_initial_call/run_PennCNV/results/res \
--out_dir ${WKDIR}/01_initial_call/run_PennCNV/results

#### (5) Merge closely adjacent CNVs and generate final results
Rscript ${WKDIR}/01_initial_call/run_PennCNV/step.5.clean.PennCNV.res.R \
--penncnv ${PENNCNV} \
--input ${WKDIR}/01_initial_call/run_PennCNV/results \
--pfb ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb

## 2 Batch effect =============================================================

### PCA on raw LRR data -------------------------------------------------------
#### (1) Randomly select 100,000 SNPs based on information from "SNP_pos.txt".
Rscript ${WKDIR}/02_batch_effect/PCA_on_LRR/step.1.down.sampling.R \
${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt \
${WKDIR}/02_batch_effect/PCA_on_LRR

#### (2) Extract LRR values at randomly selected SNPs across individuals from final report
perl ${WKDIR}/02_batch_effect/PCA_on_LRR/step.2.LRR.matrix.pl \
${WKDIR}/02_batch_effect/PCA_on_LRR/snps.down.sample.txt \
${WKDIR}/data/final_report.txt \
${WKDIR}/02_batch_effect/PCA_on_LRR/LRR_matrix_for_PCA.txt

#### (3) PCA on LRR matrix
Rscript ${WKDIR}/02_batch_effect/PCA_on_LRR/step.3.LRR.pca.R \
${WKDIR}/02_batch_effect/PCA_on_LRR \
LRR_matrix_for_PCA.txt

### PCA on summary statistics -------------------------------------------------
#### (1) Generate iPattern, PennCNV and QuantiSNP sample-level summary statistics
Rscript ${WKDIR}/02_batch_effect/PCA_on_summary_stats/step.1.prepare.stats.R \
${WKDIR}/01_initial_call/run_PennCNV/results \
${WKDIR}/02_batch_effect/PCA_on_summary_stats

#### (2) PCA on sample-level summary statistics
Rscript ${WKDIR}/02_batch_effect/PCA_on_summary_stats/step.2.stats.PCA.R \
${WKDIR}/02_batch_effect/PCA_on_summary_stats


## 3 Create CNVR ==============================================================

#### (1) Extract CNV information from individual calls made by iPattern, PennCNV and QuantiSNP
Rscript ${WKDIR}/03_create_CNVR/step.1.CNV.data.R \
${WKDIR}/03_create_CNVR \
${WKDIR}/01_initial_call/run_PennCNV/results/CNV.PennCNV_new.txt \
${WKDIR}/data/Samples_Table.txt

#### (2) Merge CNV calls from individual methods into CNVRs
Rscript ${WKDIR}/03_create_CNVR/step.2.create.CNVR.R \
--pcnv ${WKDIR}/03_create_CNVR/cnv.penncnv.txt \
--snp ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt \
--centromere ${WKDIR}/data/centromere.txt \
--output ${WKDIR}/03_create_CNVR


## 4 CNV genotyping for each CNVR =============================================

#### link necessary input data in the directory: ${WKDIR}/04_CNV_genotype/data/
cp ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb ${WKDIR}/04_CNV_genotype/data/SNP.pfb
cp ${WKDIR}/03_create_CNVR/cnvr_clean.txt ${WKDIR}/04_CNV_genotype/data/cnvr_clean.txt
cp ${WKDIR}/03_create_CNVR/cnv_clean.txt ${WKDIR}/04_CNV_genotype/data/cnv_clean.txt
cp ${WKDIR}/01_initial_call/run_PennCNV/results/CNV.PennCNV_qc_new.txt ${WKDIR}/04_CNV_genotype/data/samples_QC.txt

#### (1) Split CNVRs into different batches
Rscript ${WKDIR}/04_CNV_genotype/step.1.split.cnvrs.into.batches.R \
-i ${WKDIR}/04_CNV_genotype/data/cnvr_clean.txt \
-o ${WKDIR}/04_CNV_genotype/data/cnvr_batch.txt \
-n 200


#### (2) Submit parallelized jobs for CNV genotyping, each corresponding to one batch
Rscript ${WKDIR}/04_CNV_genotype/step.2.submit.jobs.R \
--type 0 \
--script ${WKDIR}/04_CNV_genotype \
--sourcefile ${WKDIR}/04_CNV_genotype/scripts \
--datapath ${WKDIR}/04_CNV_genotype/data \
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \
--resultpath ${WKDIR}/04_CNV_genotype/results \
--joblog ${WKDIR}/04_CNV_genotype/results

#### (3) Check submitted jobs and resubmit failed jobs
Rscript ${WKDIR}/04_CNV_genotype/step.3.check.and.resubmit.jobs.R \
--flag 1 \
--script ${WKDIR}/04_CNV_genotype \
--sourcefile ${WKDIR}/04_CNV_genotype/scripts \
--datapath ${WKDIR}/04_CNV_genotype/data \
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \
--resultpath ${WKDIR}/04_CNV_genotype/results \
--joblog ${WKDIR}/04_CNV_genotype/results

#### (4) Combine results from parallelized jobs
Rscript ${WKDIR}/04_CNV_genotype/step.4.prediction.results.R \
--datapath ${WKDIR}/04_CNV_genotype/data \
--resultpath ${WKDIR}/04_CNV_genotype/results


## 5 Boundary refinement ======================================================

#### link necessary input data in the directory: ${WKDIR}/05_boundary_refinement/data
cp ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt ${WKDIR}/05_boundary_refinement/data/SNP_pos.txt
cp ${WKDIR}/04_CNV_genotype/results/cnvr_genotype.txt ${WKDIR}/05_boundary_refinement/data/cnvr_genotype.txt
cp ${WKDIR}/04_CNV_genotype/results/matrix_CN.rds ${WKDIR}/05_boundary_refinement/data/matrix_CN.rds
cp ${WKDIR}/04_CNV_genotype/results/matrix_GQ.rds ${WKDIR}/05_boundary_refinement/data/matrix_GQ.rds

#### (1) Select CNVRs with common CNV genotype to be refined.
Rscript ${WKDIR}/05_boundary_refinement/step.1.common.CNVR.to.refine.R \
--datapath ${WKDIR}/05_boundary_refinement/data \
--resultpath ${WKDIR}/05_boundary_refinement/results \
--freq 0.05

#### (2) Submit parallelized jobs for boundary refinement, each corresponding to CNVRs in one chromosome
Rscript ${WKDIR}/05_boundary_refinement/step.2.submit.jobs.R \
--refinescript ${WKDIR}/05_boundary_refinement/CNVR.boundary.refinement.R \
--rcppfile ${WKDIR}/05_boundary_refinement/refine.cpp \
--datapath ${WKDIR}/05_boundary_refinement/data \
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \
--resultpath ${WKDIR}/05_boundary_refinement/results \
--centromere ${WKDIR}/data/centromere.txt \
--plot

#### (3) Combine results from parallelized jobs
Rscript ${WKDIR}/05_boundary_refinement/step.3.clean.results.R \
--resultpath ${WKDIR}/05_boundary_refinement/results

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## 05a_regenotype_after_refinement --------------------------------------------
## The CNVRs listed in "cnvr_regenotype_after_refine.txt" will go through CNV genotyping
## as done in 04_CNV_genotype

#### link necessary input data in the directory: ${WKDIR}/05a_regenotype_after_refinement/data/
cp ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb ${WKDIR}/05a_regenotype_after_refinement/data/SNP.pfb
cp ${WKDIR}/05_boundary_refinement/results/cnvr_regenotype_after_refine.txt ${WKDIR}/05a_regenotype_after_refinement/data/cnvr_clean.txt
cp ${WKDIR}/03_create_CNVR/cnv_clean.txt ${WKDIR}/05a_regenotype_after_refinement/data/cnv_clean.txt
cp ${WKDIR}/01_initial_call/run_PennCNV/results/CNV.PennCNV_qc_new.txt ${WKDIR}/05a_regenotype_after_refinement/data/samples_QC.txt

## the scripts for regenotyping is the same as those used in ${WKDIR}/04_CNV_genotype
#### (1) Split CNVRs into different batches
Rscript ${WKDIR}/04_CNV_genotype/step.1.split.cnvrs.into.batches.R \
-i ${WKDIR}/05a_regenotype_after_refinement/data/cnvr_clean.txt \
-o ${WKDIR}/05a_regenotype_after_refinement/data/cnvr_batch.txt \
-n 200

#### (2) Submit parallelized jobs for CNV genotyping, each corresponding to one batch
Rscript ${WKDIR}/04_CNV_genotype/step.2.submit.jobs.R \
--type 0 \
--script ${WKDIR}/04_CNV_genotype \
--sourcefile ${WKDIR}/04_CNV_genotype/scripts \
--datapath ${WKDIR}/05a_regenotype_after_refinement/data \
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \
--resultpath ${WKDIR}/05a_regenotype_after_refinement/results \
--joblog ${WKDIR}/05a_regenotype_after_refinement/results


#### (3) Check submitted jobs and resubmit failed jobs
Rscript ${WKDIR}/04_CNV_genotype/step.3.check.and.resubmit.jobs.R \
--flag 1 \
--script ${WKDIR}/04_CNV_genotype \
--sourcefile ${WKDIR}/04_CNV_genotype/scripts \
--datapath ${WKDIR}/05a_regenotype_after_refinement/data \
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \
--resultpath ${WKDIR}/05a_regenotype_after_refinement/results \
--joblog ${WKDIR}/05a_regenotype_after_refinement/results


#### (4) Combine results from parallelized jobs
Rscript ${WKDIR}/04_CNV_genotype/step.4.prediction.results.R \
--datapath ${WKDIR}/05a_regenotype_after_refinement/data \
--resultpath ${WKDIR}/05a_regenotype_after_refinement/results

#### (5) Update CN and GQ matrices as well as CNVR information
Rscript ${WKDIR}/05_boundary_refinement/step.4.update.genotype.matrix.R \
--matrixbeforerefine ${WKDIR}/05_boundary_refinement/data \
--matrixrefine ${WKDIR}/05a_regenotype_after_refinement/results \
--refinepath ${WKDIR}/05_boundary_refinement/results \
--output ${WKDIR}/05_boundary_refinement/results

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## 06_performance_assessment --------------------------------------------
mkdir -p ${WKDIR}/06_performance_assessment/results

# step1.define cnvr type
python ${WKDIR}/06_performance_assessment/step1.determine_cnvr_type.py \
${WKDIR}/05_boundary_refinement/results/cnvr_final.txt \
${WKDIR}/03_create_CNVR/cnv_create.txt \
${WKDIR}/06_performance_assessment/results/cnvr_type.txt

# step2.clean cnvr
python ${WKDIR}/06_performance_assessment/step2.clean_cnvr.py \
${WKDIR}/06_performance_assessment/results/cnvr_type.txt \
${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt

# step3.clean cnv
python ${WKDIR}/06_performance_assessment/step3.clean_cnv.py \
${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
${WKDIR}/03_create_CNVR/cnv_create.txt \
${WKDIR}/06_performance_assessment/results/final_cnv.txt \
${WKDIR}/06_performance_assessment/results/final_cnv_filted_out.txt

# step4.cnv statistics
python ${WKDIR}/06_performance_assessment/step4.cnv_stats.py \
${WKDIR}/06_performance_assessment/results/final_cnv.txt \
${WKDIR}/06_performance_assessment/results/final_cnv_statistics.csv

# step4a.initial cnv statistics
python ${WKDIR}/06_performance_assessment/step4.cnv_stats.py \
${WKDIR}/03_create_CNVR/cnv_create.txt \
${WKDIR}/06_performance_assessment/results/initial_cnv_statistics.csv

# step5.encode cnvr
python ${WKDIR}/06_performance_assessment/step5.encode_cnvr.py \
${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
${WKDIR}/06_performance_assessment/results/final_cnv.txt \
${WKDIR}/06_performance_assessment/results/inconsistent_cnvr_samples.txt \
${WKDIR}/06_performance_assessment/results/encoded_results.txt \
${WKDIR}/06_performance_assessment/results/final_samples.txt

# step6.plot cnvr distribution
Rscript ${WKDIR}/06_performance_assessment/step6.plot_cnvr_distribution.R \
-c ${PIPELINE}/01_Data/data/chromosome_size.xlsx \
-n ${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
-o ${WKDIR}/06_performance_assessment/results/cnvr_distribution.png \
-s ${WKDIR}/06_performance_assessment/results/cnvr_distribution_stats.csv

# step6a.plot cnvr frequency
Rscript ${WKDIR}/06_performance_assessment/step6a.plot_frequency_spectrum.R \
-c ${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
-e ${WKDIR}/06_performance_assessment/results/encoded_results.txt \
-d ${WKDIR}/06_performance_assessment/results/8a_cnvr_mac_data_with_id.txt \
-o ${WKDIR}/06_performance_assessment/results/8a_cnvr_mac_plot.png

# step6b.plot cnvr frequency
Rscript ${WKDIR}/06_performance_assessment/step6b.plot_frequency_spectrum_15.R \
-c ${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
-e ${WKDIR}/06_performance_assessment/results/encoded_results.txt \
-d ${WKDIR}/06_performance_assessment/results/8b_cnvr_mac_data_with_id.txt \
-o ${WKDIR}/06_performance_assessment/results/8b_cnvr_mac_plot.png

# step7.clean snp
python ${WKDIR}/06_performance_assessment/step7.clean_snp.py \
${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
${WKDIR}/data/final_report.txt \
${WKDIR}/06_performance_assessment/results/final_samples.txt \
${WKDIR}/06_performance_assessment/results/cleaned_snp.txt

# step7a.check cleaned snp
python ${WKDIR}/06_performance_assessment/step7a.check_cleaned_snp.py \
${WKDIR}/06_performance_assessment/results/cleaned_snp.txt

# step8.cnvr statistics
python ${WKDIR}/06_performance_assessment/step8.cnvr_stats.py \
${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
${WKDIR}/06_performance_assessment/results/cleaned_snp.txt \
${WKDIR}/06_performance_assessment/results/final_cnvr_statistics.csv

#****************************************PART5 Run CNV-GWAS****************************************

# Define function to run CNV-GWAS
run_cnv_gwas() {
  local trait=$1
  mkdir -p ${PIPELINE}/05_CNV_GWAS/results/${trait}
  cd ${PIPELINE}/05_CNV_GWAS/results/${trait}

  # Extract phenotraits information
  python ${PIPELINE}/05_CNV_GWAS/scripts/step1.extract_phenotypes_parityC.py \
  ${PHENOTYPE} \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/phenotypes.txt \
  ${trait} \
  ${WKDIR}/06_performance_assessment/results/final_samples.txt

  # Extract rawcnv
  bash ${PIPELINE}/05_CNV_GWAS/scripts/step2.extract_cnv_from_rawcnv.sh \
  ${WKDIR}/01_initial_call/run_PennCNV/results/CNV.PennCNV.rawcnv \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/all_formatted_cnv.txt

  # Extract cnv with phenotype
  python ${PIPELINE}/05_CNV_GWAS/scripts/step3.process_cnv_with_pheno.py \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/all_formatted_cnv.txt \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/phenotypes.txt \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/filtered_cnv.txt

  # Extract cnv with phenotrait in cnvr and have phenotype
  python ${PIPELINE}/05_CNV_GWAS/scripts/step4.process_cnv_with_pheno_in_cnvr.py \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/filtered_cnv.txt \
  ${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/filtered_cnv_in_cnvr.txt

  # Extract snp map
  python ${PIPELINE}/05_CNV_GWAS/scripts/step5.extract_snp_map.py \
  ${PIPELINE}/01_Data/results/new_snp_map_16022022_Ovine50K_4.txt \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/new_snp_map.txt

  # Run cnv-GWAS
  Rscript ${PIPELINE}/05_CNV_GWAS/scripts/step6.run_cnvgwas.R \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/filtered_cnv_in_cnvr.txt \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/phenotypes.txt \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/new_snp_map.txt \
  ${PIPELINE}/01_Data/data/chromosome_size.xlsx \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/segs_pvalue_gr_corrected.txt
}

# Run CNV-GWAS for milk, fat, and prot
for trait in milk fat prot; do
  run_cnv_gwas $trait
done

# Run phenotype statistic
python ${PIPELINE}/05_CNV_GWAS/scripts/step7.phenotype_statistic.py \
${PHENOTYPE} \
${PIPELINE}/05_CNV_GWAS/results/phenotype_statistics.txt