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
