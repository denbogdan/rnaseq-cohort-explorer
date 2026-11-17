#!/usr/bin/env python3
"""Download public TCGA RNA-seq and overall-survival data from the cBioPortal API.

For every cohort listed in config/cohorts.csv this writes two files to data/:

  <cohort>_expression.tsv.gz  genes (rows) x patients (columns); RSEM, batch-normalised
  <cohort>_survival.tsv       patient_id, os_months, death (1 = deceased, 0 = censored)

Only primary-tumour samples (TCGA sample type 01) are kept, one per patient.
Usage:  python data-raw/fetch_cohorts.py [--cohorts TCGA-LUAD] [--out data] [--max-genes N]
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://www.cbioportal.org/api"
CHUNK = 500  # genes per request


def get_json(session, method, url, retries=4, **kw):
    for attempt in range(retries):
        try:
            r = session.request(method, url, timeout=180, **kw)
            r.raise_for_status()
            return r.json()
        except requests.RequestException:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def protein_coding_genes(session):
    genes = get_json(session, "GET", f"{API}/genes",
                     params={"projection": "SUMMARY", "pageSize": 100000})
    return {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes if g.get("type") == "protein-coding"}


def fetch_expression(session, study_id, genes, max_genes=None):
    profile = f"{study_id}_rna_seq_v2_mrna"
    entrez = list(genes)[:max_genes]
    frames = []
    for i in range(0, len(entrez), CHUNK):
        body = {"entrezGeneIds": entrez[i:i + CHUNK], "sampleListId": profile}
        data = get_json(session, "POST", f"{API}/molecular-profiles/{profile}/molecular-data/fetch",
                        params={"projection": "SUMMARY"}, json=body)
        frames.append(pd.DataFrame(data)[["entrezGeneId", "sampleId", "value"]])
        print(f"  genes {min(i + CHUNK, len(entrez))}/{len(entrez)}", file=sys.stderr)
    long = pd.concat(frames, ignore_index=True)
    return long


def to_matrix(long, genes):
    """Long table -> gene x patient matrix (primary tumours, one sample per patient)."""
    long = long[long["sampleId"].str.endswith("-01")].copy()
    long["patient_id"] = long["sampleId"].str.slice(0, 12)
    long["gene"] = long["entrezGeneId"].map(genes)
    long = long.drop_duplicates(["gene", "patient_id"])
    mat = long.pivot(index="gene", columns="patient_id", values="value")
    return mat.dropna(how="all").sort_index()


def fetch_survival(session, study_id):
    rows = get_json(session, "GET", f"{API}/studies/{study_id}/clinical-data",
                    params={"clinicalDataType": "PATIENT", "projection": "SUMMARY", "pageSize": 100000})
    wide = pd.DataFrame(rows).pivot(index="patientId", columns="clinicalAttributeId", values="value")
    out = pd.DataFrame({
        "os_months": pd.to_numeric(wide["OS_MONTHS"], errors="coerce"),
        "death": wide["OS_STATUS"].str.startswith("1").astype("Int64"),
    })
    out.index.name = "patient_id"
    return out.dropna()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/cohorts.csv")
    ap.add_argument("--out", default="data")
    ap.add_argument("--cohorts", nargs="*", help="subset of cohort ids (default: all in config)")
    ap.add_argument("--max-genes", type=int, help="debug: only fetch the first N genes")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(args.config) as fh:
        cohorts = list(csv.DictReader(fh))
    if args.cohorts:
        cohorts = [c for c in cohorts if c["cohort"] in args.cohorts]

    session = requests.Session()
    genes = protein_coding_genes(session)
    for c in cohorts:
        print(f"{c['cohort']} ({c['study_id']})", file=sys.stderr)
        mat = to_matrix(fetch_expression(session, c["study_id"], genes, args.max_genes), genes)
        surv = fetch_survival(session, c["study_id"])
        keep = mat.columns.intersection(surv.index)  # patients with expression AND survival
        mat[keep].round(4).to_csv(out / f"{c['cohort']}_expression.tsv.gz", sep="\t", index_label="gene")
        surv.loc[keep].to_csv(out / f"{c['cohort']}_survival.tsv", sep="\t")
        print(f"  wrote {mat.shape[0]} genes x {len(keep)} patients "
              f"({int(surv.loc[keep, 'death'].sum())} deaths)", file=sys.stderr)


if __name__ == "__main__":
    main()
