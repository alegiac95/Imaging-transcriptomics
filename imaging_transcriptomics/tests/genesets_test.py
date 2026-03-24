from pathlib import Path

from imaging_transcriptomics.genesets import get_geneset


def test_get_geneset_resolves_packaged_and_passthrough_names():
    lake = get_geneset("lake")
    pooled = get_geneset("POOLed")
    custom = get_geneset("GO_Biological_Process_2017")

    assert Path(lake).exists()
    assert Path(pooled).exists()
    assert custom == "GO_Biological_Process_2017"
