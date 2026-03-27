from __future__ import annotations

from hashlib import sha256

import numpy as np
import pytest

import imaging_transcriptomics as imt
import imaging_transcriptomics.api as api


ATLAS_GOLDENS = {
    "dk": {
        "n_regions": 41,
        "n_genes": 15677,
        "label_hash": "d92356c388bb5e0d",
        "gene_hash": "3df2a85878a343a9",
        "value_hash": "750aecbcdce53af8",
        "corr_hash": "e6d735763d3f06a2",
    },
    "schaefer-100": {
        "n_regions": 50,
        "n_genes": 15677,
        "label_hash": "fa545bc4262ec06d",
        "gene_hash": "3df2a85878a343a9",
        "value_hash": "ec4f4968e7185463",
        "corr_hash": "acecd739e0f92a88",
    },
    "schaefer-200": {
        "n_regions": 100,
        "n_genes": 15677,
        "label_hash": "49825e5cb83e4c4b",
        "gene_hash": "3df2a85878a343a9",
        "value_hash": "db052cdc6b0b8495",
        "corr_hash": "09d4aab318a73559",
    },
    "schaefer-400": {
        "n_regions": 200,
        "n_genes": 15675,
        "label_hash": "b562fc94de9c514d",
        "gene_hash": "46cc8301fc02b910",
        "value_hash": "7418aa9fa7dec97d",
        "corr_hash": "c1179da3f1e9d67d",
    },
    "destrieux": {
        "n_regions": 74,
        "n_genes": 15675,
        "label_hash": "946b32c1e4e7a15f",
        "gene_hash": "46cc8301fc02b910",
        "value_hash": "87272b701199ac83",
        "corr_hash": "36752ff308752a5e",
    },
    "glasser-360": {
        "n_regions": 180,
        "n_genes": 15675,
        "label_hash": "22c8034debd37393",
        "gene_hash": "46cc8301fc02b910",
        "value_hash": "aefde8c9aeb6426a",
        "corr_hash": "c1179da3f1e9d67d",
    },
}

SHARED_GENE_GROUPS = {
    "genes-ahba-15677.npy": {"dk", "schaefer-100", "schaefer-200"},
    "genes-ahba-15675.npy": {"schaefer-400", "destrieux", "glasser-360"},
}


def _hash_bytes(data: bytes) -> str:
    return sha256(data).hexdigest()[:16]


def _hash_selection(selection: imt.AtlasSelection) -> tuple[str, str, str]:
    label_hash = _hash_bytes(selection.labels.to_csv(index=False, lineterminator="\n").encode())
    gene_hash = _hash_bytes("\n".join(selection.gene_labels[:, 0].astype(str)).encode())
    values = selection.expression.iloc[:, 2:].to_numpy(dtype=np.float32, copy=False)
    value_hash = _hash_bytes(values.tobytes())
    return label_hash, gene_hash, value_hash


def _hash_corr_result(result: imt.CorrelationResult) -> str:
    top = result.gene_table.head(20)
    digest = sha256()
    digest.update(top["gene"].astype(str).str.cat(sep="\n").encode())
    numeric = np.round(top[["score", "p", "fdr"]].to_numpy(dtype=np.float64), 6).astype(np.float32)
    digest.update(numeric.tobytes())
    return digest.hexdigest()[:16]


def _fake_permutations(extracted, n_permutations, *, null_method="auto", seed=1234):
    del seed
    values = extracted.values - np.mean(extracted.values)
    std = np.std(values, ddof=1)
    zvalues = values / std if std else values
    return np.tile(zvalues.reshape(-1, 1), (1, n_permutations)), null_method


def test_packaged_atlases_share_gene_label_files():
    for shared_name, atlases in SHARED_GENE_GROUPS.items():
        shared_paths = {imt.get_atlas(atlas).gene_labels_path for atlas in atlases}
        assert len(shared_paths) == 1
        shared_path = shared_paths.pop()
        assert shared_path is not None
        assert shared_path.name == shared_name
        assert shared_path.exists()


@pytest.mark.parametrize("atlas", list(ATLAS_GOLDENS), ids=list(ATLAS_GOLDENS))
def test_packaged_atlas_selection_and_corr_outputs_are_golden(monkeypatch, atlas):
    expected = ATLAS_GOLDENS[atlas]
    selection = imt.select_atlas_data(atlas=atlas, hemisphere="left", regions="all")

    assert selection.n_regions == expected["n_regions"]
    assert selection.gene_labels.shape[0] == expected["n_genes"]
    assert _hash_selection(selection) == (
        expected["label_hash"],
        expected["gene_hash"],
        expected["value_hash"],
    )

    monkeypatch.setattr(api, "_permute_scan_values", _fake_permutations)
    scan = np.linspace(-1.0, 1.0, selection.n_regions)
    result = imt.run_corr(
        scan,
        atlas=atlas,
        hemisphere="left",
        regions="all",
        n_permutations=8,
        run_gsea=False,
    )

    assert result.gene_table.shape[0] == expected["n_genes"]
    assert _hash_corr_result(result) == expected["corr_hash"]
