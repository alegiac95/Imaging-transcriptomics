import sys
from pathlib import Path

from imaging_transcriptomics.genesets import (
    get_geneset,
    list_packaged_genesets,
    list_remote_genesets,
    resolve_geneset_resource,
)


def test_get_geneset_resolves_packaged_and_passthrough_names():
    lake = get_geneset("lake")
    pooled = get_geneset("POOLed")
    custom = get_geneset("GO_Biological_Process_2017")

    assert Path(lake).exists()
    assert Path(pooled).exists()
    assert custom == "GO_Biological_Process_2017"


def test_list_packaged_genesets_contains_expected_names():
    assert list_packaged_genesets() == ("lake", "pooled")


def test_resolve_geneset_resource_downloads_remote_library(monkeypatch):
    fake_gseapy = type(
        "FakeGSEApy",
        (),
        {
            "get_library": staticmethod(
                lambda name, organism="Human", min_size=0, max_size=100_000: {"Term": ["GENE1", "GENE2"]}
            ),
            "get_library_name": staticmethod(lambda organism="Human": ["RemoteA", "RemoteB"]),
        },
    )
    monkeypatch.setitem(sys.modules, "gseapy", fake_gseapy)

    resolved = resolve_geneset_resource("RemoteA", organism="Mouse")

    assert resolved == {"Term": ["GENE1", "GENE2"]}
    assert list_remote_genesets("Mouse") == ["RemoteA", "RemoteB"]
