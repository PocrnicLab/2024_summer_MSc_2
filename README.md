# Copy number variant genome-wide association study for milk production traits in sheep
This work was completed as part of a master's thesis for a degree in Bioinformatics at the University of Edinburgh, conducted from May to August 2024. The thesis was undertaken by B245390 at the Roslin Institute under the guidance of Dr. Ivan Pocrnic and Dr. Chrissy Rochus.

# Purposes
This bioinformatics pipeline is designed to analyze SNP data from sheep genotyped with the Illumina OvineSNP50 BeadChip. The pipeline includes data quality control using PLINK, CNV identification with PennCNV, and CNV region (CNVR) analysis via EnsembleCNV. Genome-wide association studies (GWAS) on milk production traits are conducted using both SNPs and CNVs. Functional annotation is performed to identify candidate genes and quantitative trait loci (QTLs) associated with daily milk, fat, and protein production.
![RP_flowchart_final (1)](https://github.com/user-attachments/assets/6803b86c-e8fc-4d76-ab19-7bdb8184ed02)

# Computational environment setting
The development and testing of the pipeline were conducted on the University of Edinburgh's compute cluster, Eddie Mark 3, which utilizes the Altair Gridengine batch system on Scientific Linux 7. The software packages R, Python, and PLINK were loaded on the Eddie system using the following module commands: 
‘module load python/3.11.4’
‘module load R/4.3.0’
‘module load roslin/plink/1.90p’.

PennCNV.1.0.5 was downloaded and installed following the instructions at 
https://penncnv.openbioinformatics.org/en/latest/user-guide/install/. 

EnsembleCNV was installed from 
https://github.com/HaoKeLab/ensembleCNV, 
with code modifications made according to the CNV identification software used and operating environment. Additionally, the versions of third-party R packages and Python libraries used are listed below:

Package Name	Programming Language	Version

CNVRanger	R	1.18.1

cowplot	R	1.1.3

dplyr	R	1.1.4

GALLO	R	1.1

ggplot2	R	3.5.1

gridExtra	R	2.3

matplotlib	Python	3.9.0

mclust	R	6.0.1

mixtools	R	2.0.0

modeest	R	2.4.0

numpy	Python	2.0.0

optparse	R	1.7.5

pandas	Python	2.0.3

pheatmap	R	1.0.12

plyr	R	1.8.9

qqman	R	0.1.9

RColorBrewer	R	1.1.3

Rcpp	R	1.0.12

readxl	R	1.4.3

requests	Python	2.32.3

scales	R	1.3.0

scipy	Python	1.14.0

statsmodels	Python	0.14.2

tibble	R	3.2.1

tidyr	R	1.3.1

tools	R	4.3.2

The following part of https://github.com/PocrnicLab/2024_summer_MSc_2/blob/main/main.sh sets up the library and installs necessary packages, which can be deleted if these steps are already done.
#=====================================================================================

cd ${PIPELINE}

module load python/3.11.4
module load R/4.3.0

R_LIBS_DIR="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/R_libs/"

export R_LIBS=$R_LIBS_DIR

Rscript -e "
packages <- c('cowplot', 'data.table', 'dplyr', 'ggplot2', 'gridExtra', 'mclust', 'mixtools', 'modeest', 'optparse', 'pheatmap', 'plyr', 'RColorBrewer', 'Rcpp', 'RcppArmadillo', 'tibble')
install_if_missing

PY_LIB_PATH="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/py_libs/"
export PYTHONPATH="${PY_LIB_PATH}:$PYTHONPATH"

REQUIRED_PKG=("pandas" "matplotlib" "seaborn" "scipy" "numpy" "statsmodels")

install_packages() 

install_packages
#=====================================================================================

# Input files preparation
The preparation of input files is essential for the operation of the pipeline. In addition to the Illumina final report, the corresponding SNP map and phenotype records, the following files are also required: 
1. centromere.txt, which includes the columns  ‘chr’ for chromosome and ‘position’, along with the relevant data. This file is used in EnsembleCNV to determine whether a CNVR is located on the p-arm or q-arm of the chromosome. Due to the difficulty in obtaining centromere data for sheep, I used a “pseudo” centromere file in this experiment, setting the position of all centromeres to 0 bp. Consequently, all CNVRs will be classified as being on the q-arm.
2. chromosome_size.xlsx, which must at least include the columns ‘Chromosome’ and ‘Size (bp)’ along with the corresponding data. This file is used to plot the CNVR distribution by chromosomes. The relevant data can be obtained from NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/) and saved as a spreadsheet.
3.  ‘.gtf’ files for gene annotation, which can be found in ensembl FTP website (https://www.ensembl.org/info/data/ftp/index.html).
4.  ‘.gff’ files for quantitative trait loci (QTLs) annotation, which can be found in Animal QTLdb (https://www.animalgenome.org/QTLdb/).

It is noteworthy that all input files should use the same genome assembly version. While the pipeline can still run with input files of different assembly versions, the results may not be biologically sensible. In this experiment, I used the easily assessable Genome assembly Oar_v3.1 data. Since the Illumina final report and the corresponding SNP map were based on different assembly versions, an additional mapping file was downloaded from 
https://webserver.ibba.cnr.it/SNPchimp/ 
to update these two files. In the mapping file, SNPs deleted in the target assembly version are recorded as being located at the 0 bp position on chromosome 99.

All external data files except for '.gtf' file are available at https://github.com/PocrnicLab/2024_summer_MSc_2/tree/main/01_Data/data

# How to use this pipeline
After define the needed pathway and threshold in main.sh, for example:

#=====================================================================================
PIPELINE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2"
PENNCNV="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/02_Tools_and_Computational_Environment/PennCNV-1.0.5/"
FINALREPORT="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_final_report_16022022_Ovine50K_4.txt"
SNPMAP="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/we_uz_snp_map_16022022_Ovine50K_4.txt"
PHENOTYPE="/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/01_Data/Pag_sheep_phenotypes_OneRecordPerParity.txt"
PVALUEFORCNV=0.10
PVALUEFORSNP=0.05
#=====================================================================================

The whole pipeline can be run easily by typing:

bash main.sh
