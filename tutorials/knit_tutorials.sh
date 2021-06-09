#!/bin/bash
HL1='\033[0;32m'
NC='\033[0m' # No Color

printf "\n${HL1} (1) -----------------------------------------------------------------\n"
date
printf "${NC}\n"
Rscript -e 'library(rmarkdown); rmarkdown::render("basic/basic.Rmd", "html_document")'

printf "\n${HL1} (2) -----------------------------------------------------------------\n"
date
printf "${NC}\n"
Rscript -e 'library(rmarkdown); rmarkdown::render("heter/heter.Rmd", "html_document")'

printf "\n${HL1} (3) -----------------------------------------------------------------\n"
date
printf "${NC}\n"
Rscript -e 'library(rmarkdown); rmarkdown::render("uncrt/uncrt.Rmd", "html_document")'

printf "\n${HL1} (4) -----------------------------------------------------------------\n"
date
printf "${NC}\n"
Rscript -e 'library(rmarkdown); rmarkdown::render("bb/bb.Rmd", "html_document")'
date

printf "\n${HL1} (FINISH) ------------------------------------------------------------\n"
date
printf "${NC}\n"
