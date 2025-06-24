To repeat the data from this article, follow these steps:

1. Install Apptainer:
```
https://apptainer.org/docs/admin/main/installation.html
```
2. Download Apptainer container via browser:
```
https://files.sberdisk.ru/s/w8Cx8CNLzCWScgW
```
3. Run Apptainer container:
```
apptainer run r-niknit_V3.sif
```
4. Download code and data from Github:
```
wget https://github.com/niknit96/Aleynova_et.al.2025/archive/master.zip
unzip master.zip
cd ./Aleynova_et.al.2025-main
```
5. Run R scripts:
```
Rscript -e 'rmarkdown::render("./16S/16S_new.Rmd", output_file="../16S_NGS_report.html")'
Rscript -e 'rmarkdown::render("./ITS/ITS_new.Rmd", output_file="../ITS_NGS_report.html")'
```
6. Results in ./Aleynova_et.al.2025-main.
