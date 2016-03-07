all: small braya

small-data: ./analysis/small-water/small-water-age-depth-model.R
	Rscript "./analysis/small-water/small-water-age-depth-model.R"

small-example: ./analysis/small-water/small-water-analysis.R
	Rscript "./analysis/small-water/small-water-analysis.R"

small: small-data small-example

braya: ./analysis/braya-so/braya-so-analysis.R
	Rscript ./analysis/braya-so/braya-so-analysis.R
