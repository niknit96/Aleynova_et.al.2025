To repeat the data from this article, follow these steps:

1. Install Apptainer on Linux:
https://apptainer.org/docs/admin/main/installation.html
2. Download code and data from Github:
```
wget https://github.com/niknit96/Aleynova_et.al.2024/archive/master.zip
unzip master.zip
cd ./Aleynova_et.al.2025-main
```
3. Download Apptainer container:
```
wget https://biosoil.ru/downloads/biotech/Containers%20(Apptainer)/r-niknit_V3.sif
```
4. Run container and R scripts:
```
apptainer exec r-niknit_V3.sif \
  Rscript -e 'rmarkdown::render("16S.Rmd", output_file="16S_NGS_report.html")' \
  Rscript -e 'rmarkdown::render("ITS.Rmd", output_file="ITS_NGS_report.html")'
```
5. Results in ./Aleynova_et.al.2025-main.
