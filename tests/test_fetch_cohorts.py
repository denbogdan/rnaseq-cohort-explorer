import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "data-raw"))
from fetch_cohorts import to_matrix  # noqa: E402


def test_to_matrix_keeps_primary_tumours_one_per_patient():
    genes = {1: "A1BG", 2: "A2M"}
    long = pd.DataFrame({
        "entrezGeneId": [1, 2, 1, 1, 2],
        "sampleId": ["TCGA-AA-0001-01", "TCGA-AA-0001-01", "TCGA-AA-0001-11", "TCGA-AA-0002-01", "TCGA-AA-0002-01"],
        "value": [1.0, 2.0, 99.0, 3.0, 4.0],
    })
    mat = to_matrix(long, genes)
    assert list(mat.columns) == ["TCGA-AA-0001", "TCGA-AA-0002"]  # -11 (normal) dropped
    assert mat.loc["A1BG", "TCGA-AA-0001"] == 1.0
    assert mat.shape == (2, 2)
