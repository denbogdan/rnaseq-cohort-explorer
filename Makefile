.PHONY: data run test docker
data:            ## download public TCGA data (~5 min)
	python data-raw/fetch_cohorts.py
run:             ## launch the Shiny dashboard on http://localhost:3838
	Rscript -e "rmarkdown::run('app/app.Rmd', shiny_args = list(port = 3838))"
test:
	pytest tests && Rscript -e "testthat::test_dir('tests/testthat')"
docker:
	docker build -t rnaseq-cohort-explorer .
