To repeat the data from this article, follow these steps:

1. Install Apptainer on Linux:
```
https://apptainer.org/docs/admin/main/installation.html
```
3. Download Apptainer container via browser:
```
https://files.sberdisk.ru/s/w8Cx8CNLzCWScgW
```
4. Run Apptainer container:
```
apptainer run r-niknit_V3.sif
```
5. Download code and data from Github:
```
wget https://github.com/niknit96/Aleynova_et.al.2024/archive/master.zip
unzip master.zip
cd ./Aleynova_et.al.2025-main
```
4. Run R scripts:
```
Rscript -e 'rmarkdown::render("./16S/16S_new.Rmd", output_file="16S_NGS_report.html")' \
Rscript -e 'rmarkdown::render("./ITS/ITS_new.Rmd", output_file="ITS_NGS_report.html")'
```
5. Results in ./Aleynova_et.al.2025-main.
