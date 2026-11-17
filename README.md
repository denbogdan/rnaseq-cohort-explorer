# RNA-seq Cohort Explorer

An interactive R Shiny dashboard for exploring gene expression in patient RNA-seq cohorts.
Originally built for a cancer-biomarker research team to let wet-lab colleagues query internal
patient datasets without writing code; this version runs on **public TCGA data** (lung
adenocarcinoma and kidney clear cell carcinoma).

<!-- TODO: add a screenshot, e.g. docs/screenshot.png, after running the app -->

## Features

| Tab | What it does |
|---|---|
| Gene Expression | Spearman correlation of two genes across patients, interactive scatter plot, per-patient expression bars |
| Pathway Analysis | Pre-ranked GSEA (`clusterProfiler`) of genes ranked by correlation with a query gene, against MSigDB Hallmark / Reactome / KEGG; downloadable table; results cached on disk |
| Overall Survival | Kaplan-Meier curves split by expression, using the optimal max-stat cutpoint or the median |

## Quick start

```bash
git clone <this repo> && cd rnaseq-cohort-explorer
pip install -r requirements.txt
python data-raw/fetch_cohorts.py       # downloads public TCGA data from cBioPortal into data/
Rscript install_packages.R
make run                               # http://localhost:3838
```

The first launch builds the MSigDB gene-set cache (`data/genesets.rds`) with `msigdbr`.

### Docker

```bash
python data-raw/fetch_cohorts.py
docker build -t rnaseq-cohort-explorer .
docker run -p 3838:3838 -v $(pwd)/data:/srv/app/data rnaseq-cohort-explorer
```

## Layout

```
app/app.Rmd            flexdashboard + Shiny UI (thin layer)
R/analysis.R           pure analysis functions: correlation, ranking, GSEA, survival stratification
R/data.R               data and gene-set loading
data-raw/              reproducible download of public data
config/cohorts.csv     cohorts to include (add a row to add a cohort)
tests/                 testthat (R) and pytest (Python) unit tests
```

Adding a cohort: append a `cohort,label,study_id` row to `config/cohorts.csv` (any cBioPortal
study with an `rna_seq_v2_mrna` profile and `OS_MONTHS`/`OS_STATUS` clinical data), then re-run
the fetch script.

## Data and licence

Data: TCGA PanCancer Atlas via cBioPortal (public). Please follow TCGA/cBioPortal citation guidance.
Code: MIT.

## Notes on the original

Differences from the internal version: private cohorts replaced by public ones; analysis logic
moved out of the Rmd into testable functions; cohort sizes and plot layouts no longer hard-coded;
gene symbols validated with clear errors; data fetched reproducibly instead of shipped as `.RData`;
usage logging and internal hostnames/contacts removed.
