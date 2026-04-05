from __future__ import annotations

from imaging_transcriptomics.cli import _resolve_run_gsea, build_parser
from imaging_transcriptomics.cli_support.runners import resolve_enrichment_method


class _Args:
    def __init__(self, run_gsea=None, ora_p_threshold=None, enrichment=None):
        self.run_gsea = run_gsea
        self.ora_p_threshold = ora_p_threshold
        self.enrichment = enrichment


def test_resolve_enrichment_defaults_to_ensemble():
    assert resolve_enrichment_method(_Args(run_gsea=None, ora_p_threshold=None, enrichment=None)) == "ensemble"


def test_resolve_enrichment_defaults_to_ora_when_threshold_is_given():
    assert resolve_enrichment_method(_Args(run_gsea=None, ora_p_threshold=0.05, enrichment=None)) == "ora"


def test_resolve_enrichment_respects_explicit_gsea_request():
    assert resolve_enrichment_method(_Args(run_gsea=True, ora_p_threshold=0.05, enrichment=None)) == "gsea"


def test_resolve_enrichment_respects_explicit_none_request():
    assert resolve_enrichment_method(_Args(run_gsea=False, ora_p_threshold=None, enrichment="none")) == "none"


def test_resolve_run_gsea_helper_only_reports_true_for_gsea():
    assert _resolve_run_gsea(_Args(run_gsea=None, ora_p_threshold=None, enrichment=None)) is False
    assert _resolve_run_gsea(_Args(run_gsea=True, ora_p_threshold=None, enrichment=None)) is True


def test_top_level_help_is_descriptive():
    help_text = build_parser().format_help()
    assert "built-in atlas handling" in help_text
    assert "Examples:" in help_text
    assert "imt corr" in help_text
    assert "imt gene" in help_text
    assert "imt genesets" in help_text
    assert "gene-pca" in help_text
    assert "gedar" in help_text


def test_corr_help_mentions_enrichment_choices():
    parser = build_parser()
    corr_parser = parser._subparsers._group_actions[0].choices["corr"]
    help_text = corr_parser.format_help()
    assert "run one enrichment backend" in help_text
    assert "--enrichment" in help_text
    assert "ensemble" in help_text
    assert "ora-p-threshold" in help_text
    assert "Legacy compatibility flag equivalent" in help_text
    assert "--no-gsea" in help_text
    assert "--geneset-organism" in help_text
    assert "--regions" in help_text
    assert "default" in help_text
    assert "cort+sub" in help_text


def test_gene_pca_help_mentions_gene_list_and_normalization():
    parser = build_parser()
    gene_pca_parser = parser._subparsers._group_actions[0].choices["gene-pca"]
    help_text = gene_pca_parser.format_help()
    assert "normalize those genes across regions" in help_text
    assert "--genes" in help_text
    assert "--ncomp" in help_text


def test_gene_help_mentions_expression_and_coexpression_controls():
    parser = build_parser()
    gene_parser = parser._subparsers._group_actions[0].choices["gene"]
    help_text = gene_parser.format_help()
    assert "co-expressed genes" in help_text
    assert "--gene" in help_text
    assert "--top-n" in help_text
    assert "--fdr-threshold" in help_text
    assert "--raw-expression" in help_text


def test_gedar_help_mentions_weighted_projection():
    parser = build_parser()
    gedar_parser = parser._subparsers._group_actions[0].choices["gedar"]
    help_text = gedar_parser.format_help()
    assert "weighted gene table" in help_text
    assert "--weights" in help_text
    assert "--weight-column" in help_text
    assert "--normalize-expression" in help_text
    assert "--enrichment" in help_text
    assert "--geneset" in help_text
    assert "split up/down ORA" in help_text


def test_genesets_help_mentions_packaged_and_enrichr_sources():
    parser = build_parser()
    geneset_parser = parser._subparsers._group_actions[0].choices["genesets"]
    help_text = geneset_parser.format_help()
    assert "Packaged entries are always available" in help_text
    assert "--packaged-only" in help_text
    assert "--organism" in help_text
