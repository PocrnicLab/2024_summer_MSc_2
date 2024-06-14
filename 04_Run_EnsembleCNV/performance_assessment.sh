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
#====================================================================================================
## 06 performance assessment --------------------------------------------
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
-o ${WKDIR}/06_performance_assessment/results/cnvr_distribution.png

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