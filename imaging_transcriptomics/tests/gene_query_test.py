from pathlib import Path

from imaging_transcriptomics import load_expression_frame, run_gene


def _example_gene() -> str:
    frame = load_expression_frame("dk", hemisphere="left", regions="default", zscore_expression=True)
    return str(frame.columns[2])


def test_run_gene_returns_expression_and_coexpression_tables():
    gene = _example_gene()
    result = run_gene(gene, atlas="dk", top_n=12)

    assert result.gene == gene
    assert "expression_z" in result.regional_values.columns
    assert result.regional_values.shape[0] > 0
    assert result.gene not in set(result.gene_table["gene"].astype(str))
    assert result.coexpression_matrix.index.tolist() == result.coexpression_matrix.columns.tolist()
    assert result.coexpression_matrix.index[0] == gene
    assert result.coexpressed_genes.shape[0] <= 12
    assert result.anticorrelated_genes.shape[0] <= 12
    if not result.coexpressed_genes.empty:
        assert result.coexpressed_genes["selected"].astype(bool).all()
        assert (result.coexpressed_genes["score"] > 0).all()
        assert (result.coexpressed_genes["fdr"] <= result.fdr_threshold).all()
    if not result.anticorrelated_genes.empty:
        assert result.anticorrelated_genes["selected_negative"].astype(bool).all()
        assert (result.anticorrelated_genes["score"] < 0).all()
        assert (result.anticorrelated_genes["fdr"] <= result.fdr_threshold).all()


def test_run_gene_writes_standard_bundle(tmp_path: Path):
    gene = _example_gene()
    result = run_gene(gene, atlas="dk", output_dir=tmp_path)

    assert result.output_dir == tmp_path
    assert (tmp_path / "README.txt").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "gene_query_summary.tsv").exists()
    assert (tmp_path / "gene_expression.tsv").exists()
    assert (tmp_path / "top_expression_regions.tsv").exists()
    assert (tmp_path / "top_coexpressed_genes.tsv").exists()
    assert (tmp_path / "top_negatively_correlated_genes.tsv").exists()
    assert (tmp_path / "plots" / "gene_expression_brain.png").exists()
    assert (tmp_path / "plots" / "gene_expression_cortex.png").exists()
    assert (tmp_path / "plots" / "gene_coexpression_matrix.png").exists()
    assert (tmp_path / "plots" / "gene_coexpression_top_genes.png").exists()
    assert not (tmp_path / "gene_coexpression.tsv").exists()
    assert not (tmp_path / "plots" / "gene_expression.png").exists()
    assert not (tmp_path / "plots" / "gene_coexpression_distribution.png").exists()
