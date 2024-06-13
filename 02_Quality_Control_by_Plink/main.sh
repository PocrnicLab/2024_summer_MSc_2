#!/bin/bash

PIPELINE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2"

FINALREPORT="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_final_report_16022022_Ovine50K_4.txt"
SNPMAP="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_snp_map_16022022_Ovine50K_4.txt"

#====================================================================================================
module load python/3.11.4
#====================================================================================================

mkdir -p ${PIPELINE}/02_Quality_Control_by_Plink/results

# generate files for plink
python ${PIPELINE}/02_Quality_Control_by_Plink/scripts/step1.generate_bed_map.py \
  --genotype_file "${PIPELINE}/01_Data/results/new_final_report_16022022_Ovine50K_4.txt" \
  --snp_map_file "${PIPELINE}/01_Data/results/new_snp_map_16022022_Ovine50K_4.txt" \
  --ped_file "${PIPELINE}/02_Quality_Control_by_Plink/results/input_data.ped" \
  --map_file "${PIPELINE}/02_Quality_Control_by_Plink/results/input_data.map"

#====================================================================================================
cd ${PIPELINE}/02_Quality_Control_by_Plink/results

bash ${PIPELINE}/02_Quality_Control_by_Plink/scripts/step2.run_plink.sh \
${PIPELINE}/02_Quality_Control_by_Plink/results