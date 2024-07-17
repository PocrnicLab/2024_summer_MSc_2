#!/bin/bash

## Following instructions to install ensembleCNV, PennCNV.
## Assume the installation folder
PIPELINE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2"
PENNCNV="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/PennCNV-1.0.5/"
FINALREPORT="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_final_report_16022022_Ovine50K_4.txt"
SNPMAP="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_snp_map_16022022_Ovine50K_4.txt"
PHENOTYPE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/Pag_sheep_phenotypes_OneRecordPerParity.txt"
PVALUEFORCNV=0.10
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
if (!requireNamespace('devtools', quietly = TRUE)) {
  install.packages('devtools', lib = '$R_LIBS_DIR', repos = 'http://cran.us.r-project.org')
}
library(devtools)
if (!require('GALLO', character.only = TRUE, quietly = TRUE)) {
  devtools::install_github('pablobio/GALLO')
}
"
# Set Python library folder
PY_LIB_PATH="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/py_libs/"
export PYTHONPATH="${PY_LIB_PATH}:$PYTHONPATH"

# Required packages
REQUIRED_PKG=("pandas" "matplotlib" "seaborn" "scipy" "numpy" "statsmodels")

# Function to check and install required packages
install_packages() {
    for pkg in "${REQUIRED_PKG[@]}"; do
        if ! python3 -c "import ${pkg}" &> /dev/null; then
            echo "Installing ${pkg}..."
            pip install --target="${PY_LIB_PATH}" "${pkg}"
        else
            echo "${pkg} is already installed."
        fi
    done
}

# Install required packages
install_packages
#****************************************PART6 Functional Annotation for CNV GWAS****************************************

# Define function to run functional annotation for a given trait
run_functional_annotation() {
  local trait=$1
  mkdir -p ${PIPELINE}/06_Functional_Annotation/results/${trait}
  cd ${PIPELINE}/06_Functional_Annotation/results/${trait}

  # Extract CNVR after GWAS
  python ${PIPELINE}/06_Functional_Annotation/scripts/step1.extract_cnvr_after_gwas.py \
  ${PIPELINE}/05_CNV_GWAS/results/${trait}/segs_pvalue_gr_corrected.txt \
  ${WKDIR}/06_performance_assessment/results/final_cnvr_type.txt \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/cnvr_gwas.txt \
  --pvalue_threshold ${PVALUEFORCNV}

  # Run Gallo
  Rscript ${PIPELINE}/06_Functional_Annotation/scripts/step2.run_gallo.R \
  --cnvr ${PIPELINE}/06_Functional_Annotation/results/${trait}/cnvr_gwas.txt \
  --qtl ${PIPELINE}/01_Data/data/merged_QTLdb_sheepOAR3.gff \
  --gene ${PIPELINE}/01_Data/data/Ovis_aries.Oar_v3.1.112.gtf \
  --output_dir ${PIPELINE}/06_Functional_Annotation/results/${trait}
  
  # Enrichment analysis based on QTL Name
  python ${PIPELINE}/06_Functional_Annotation/scripts/step3.enrichment_analysis.py \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/output_qtls.txt \
  ${PIPELINE}/01_Data/data/merged_QTLdb_sheepOAR3.gff \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/name_enrichment_analysis_bubble_plot.png

  # Enrichment analysis based on QTL Name - Chr
  python ${PIPELINE}/06_Functional_Annotation/scripts/step3a.enrichment_analysis_chr.py \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/output_qtls.txt \
  ${PIPELINE}/01_Data/data/merged_QTLdb_sheepOAR3.gff \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/name_chr_enrichment_analysis_bubble_plot.png
  
  # Fetch GO and KEGG term
  python ${PIPELINE}/06_Functional_Annotation/scripts/step4.fetch_go_kegg.py \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/output_genes.txt \
  ${PIPELINE}/06_Functional_Annotation/results/${trait}/go_kegg_info.csv
}

# Run functional annotation for milk, fat, and prot
for trait in milk fat prot; do
  run_functional_annotation $trait
done

