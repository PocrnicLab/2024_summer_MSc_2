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
#====================================================================================================
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
