To repeat the data from this article, follow these steps:

1. Install mamba or conda on Linux:
* mamba: https://github.com/mamba-org/mamba (more faster than conda)
* conda: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html
2. Download code and data from Github:
```
wget https://github.com/niknit96/Aleynova_et.al.2025/archive/master.zip
unzip master.zip
```
3. Create and activate conda environment for download SRA files and pre-trained classifiers for QIIME 2 Scikit-learn algorithm:
```
mamba env create --file ./Aleynova_et.al.2025-main/sra-tools.yml || conda env create --file ./Aleynova_et.al.2025-main/sra-tools.yml
mamba activate sra-tools || conda activate sra-tools
```
4. Run bash script for download SRA files and pre-trained classifiers for QIIME 2 Scikit-learn algorithm:
```
bash ./Aleynova_et.al.2025-main/download.sh
```
5. Create and activate conda environment for analysis data:
```
mamba env create --file ./Aleynova_et.al.2025-main/Aleynova_et.al.2025.yml || conda env create --file ./Aleynova_et.al.2025-main/Aleynova_et.al.2025.yml
mamba activate Aleynova_et.al.2025 || conda activate Aleynova_et.al.2025
```
6. Download and install qiime2R and microViz packages (accessed on 11 February 2025):
```
Rscript -e 'devtools::install_github("jbisanz/qiime2R")'
Rscript -e 'devtools::install_github("david-barnett/microViz")'
```
7. Run bash script for analysis data:
```
bash ./Aleynova_et.al.2025-main/analysis.sh
```
8. Results in ./Aleynova_et.al.2025-main.
