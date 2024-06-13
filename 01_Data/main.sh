#!/bin/bash

PIPELINE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/"

FINALREPORT="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_final_report_16022022_Ovine50K_4.txt"
SNPMAP="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_snp_map_16022022_Ovine50K_4.txt"

#====================================================================================================
module load python/3.11.4
#====================================================================================================

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