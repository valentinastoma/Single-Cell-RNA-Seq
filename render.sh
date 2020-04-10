#!/usr/bin/env bash
# render.sh
R -e "rmarkdown::render('SingleCellSeurat.Rmd', output_format='all')"