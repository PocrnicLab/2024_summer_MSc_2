#!/bin/bash

## Following instructions to install ensembleCNV, PennCNV.
## Assume the installation folder
PIPELINE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2"
PENNCNV="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/PennCNV-1.0.5/"
FINALREPORT="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_final_report_16022022_Ovine50K_4.txt"
SNPMAP="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_snp_map_16022022_Ovine50K_4.txt"

#====================================================================================================
ENSEMBLECNV=${PIPELINE}/04_Run_EnsembleCNV
## Working directory for a ensemble new project
WKDIR=${PIPELINE}/04a_EnsembleCNV_Working_Directory
#====================================================================================================
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

## Create a new project
cd $ENSEMBLECNV
chmod +x create_new_project.sh

./create_new_project.sh $WKDIR

mkdir -p ${PIPELINE}/03_Generate_Input_Files/results
## 1 Initial call =============================================================

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
## 6 Performance assessment ===================================================

