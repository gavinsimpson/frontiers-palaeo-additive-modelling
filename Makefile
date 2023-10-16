all: small braya correlfun

small-data: ./analysis/small-water/small-water-age-depth-model.R
	Rscript "./analysis/small-water/small-water-age-depth-model.R"

small-example: ./analysis/small-water/small-water-analysis.R
	Rscript "./analysis/small-water/small-water-analysis.R"

small: small-data small-example

braya-example: ./analysis/braya-so/braya-so-analysis.R
	Rscript ./analysis/braya-so/braya-so-analysis.R

braya-fail: ./analysis/braya-so/braya-so-failure-analysis.R
	Rscript ./analysis/braya-so/braya-so-failure-analysis.R

braya: braya-example braya-fail

correlfun: ./analysis/other/correlation-functions.R
	Rscript ./analysis/other/correlation-functions.R

manuscript:  header.tex manuscript.Rmd
	Rscript -e "library(rmarkdown); render('manuscript.Rmd')"

markdown:
	Rscript -e "library(rmarkdown); render('manuscript.Rmd', output_format = 'md_document')"

odt:
	Rscript -e "library(rmarkdown); render('manuscript.Rmd', output_format = 'odt_document')"

word:
	Rscript -e "library(rmarkdown); render('manuscript.Rmd', output_format = 'word_document')"

script: manuscript.Rmd
	Rscript -e "knitr::purl('manuscript.Rmd')"

supplement: header.tex supplementary-materials.Rmd
	Rscript -e "library(rmarkdown); render('supplementary-materials.Rmd')"

preprint: manuscript.pdf supplementary-materials.pdf
	pdfjam manuscript.pdf supplementary-materials.pdf -o preprint.pdf

