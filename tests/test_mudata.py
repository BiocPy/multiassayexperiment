import numpy as np
import pandas as pd
import pytest
import biocframe
from summarizedexperiment import SummarizedExperiment
from singlecellexperiment import SingleCellExperiment
from multiassayexperiment import MultiAssayExperiment
import mudata as mu

def test_MAE_mudata_roundtrip():
    # Construct a SummarizedExperiment
    se = SummarizedExperiment(
        assays={"counts": np.random.rand(10, 3)},
        column_data=pd.DataFrame({"group": ["A", "B", "A"]}, index=["cell_1", "cell_2", "cell_3"]),
        row_data=pd.DataFrame({"gene_symbol": [f"G_{i}" for i in range(10)]}, index=[f"gene_{i}" for i in range(10)])
    )

    # Construct a SingleCellExperiment
    sce = SingleCellExperiment(
        assays={"counts": np.random.rand(15, 2)},
        column_data=pd.DataFrame({"batch": ["1", "2"]}, index=["sc_1", "sc_2"]),
        row_data=pd.DataFrame({"gene_type": [f"T_{i}" for i in range(15)]}, index=[f"scgene_{i}" for i in range(15)])
    )

    # Construct MultiAssayExperiment
    column_data = biocframe.BiocFrame(
        {"patient": ["P1", "P2", "P3", "P4", "P5"]},
        row_names=["patient_1", "patient_2", "patient_3", "patient_4", "patient_5"]
    )
    sample_map = biocframe.BiocFrame({
        "colname": ["cell_1", "cell_2", "cell_3", "sc_1", "sc_2"],
        "assay": ["se", "se", "se", "sce", "sce"],
        "primary": ["patient_1", "patient_2", "patient_3", "patient_4", "patient_5"]
    })

    mae = MultiAssayExperiment(
        experiments={"se": se, "sce": sce},
        column_data=column_data,
        sample_map=sample_map,
        metadata={"study_name": "Multiomics Study"}
    )

    # Transform to MuData
    mdata = mae.to_mudata()
    assert isinstance(mdata, mu.MuData)
    assert "se" in mdata.mod
    assert "sce" in mdata.mod
    assert mdata.uns["study_name"] == "Multiomics Study"
    assert "sample_map" in mdata.uns

    # Verify column_data in uns
    col_df = mdata.uns["column_data"]
    assert list(col_df.index) == ["patient_1", "patient_2", "patient_3", "patient_4", "patient_5"]
    assert list(col_df["patient"]) == ["P1", "P2", "P3", "P4", "P5"]

    # Round trip from MuData
    mae_rt = MultiAssayExperiment.from_mudata(mdata)
    assert isinstance(mae_rt, MultiAssayExperiment)
    assert "se" in mae_rt.experiments
    assert "sce" in mae_rt.experiments
    assert list(mae_rt.column_data.row_names) == ["patient_1", "patient_2", "patient_3", "patient_4", "patient_5"]
    assert list(mae_rt.sample_map["colname"]) == ["cell_1", "cell_2", "cell_3", "sc_1", "sc_2"]
    assert mae_rt.metadata["study_name"] == "Multiomics Study"
