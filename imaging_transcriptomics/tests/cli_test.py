from __future__ import annotations

from imaging_transcriptomics.script.imagingtranscriptomics import _resolve_run_gsea, build_parser


class _Args:
    def __init__(self, run_gsea=None, ora_p_threshold=None):
        self.run_gsea = run_gsea
        self.ora_p_threshold = ora_p_threshold


def test_resolve_run_gsea_defaults_to_true_without_ora():
    assert _resolve_run_gsea(_Args(run_gsea=None, ora_p_threshold=None)) is True


def test_resolve_run_gsea_defaults_to_false_for_ora_only_runs():
    assert _resolve_run_gsea(_Args(run_gsea=None, ora_p_threshold=0.05)) is False


def test_resolve_run_gsea_respects_explicit_enable_with_ora():
    assert _resolve_run_gsea(_Args(run_gsea=True, ora_p_threshold=0.05)) is True


def test_resolve_run_gsea_respects_explicit_disable():
    assert _resolve_run_gsea(_Args(run_gsea=False, ora_p_threshold=None)) is False


def test_top_level_help_is_descriptive():
    help_text = build_parser().format_help()
    assert "built-in atlas handling" in help_text
    assert "Examples:" in help_text
    assert "imagingtranscriptomics corr" in help_text


def test_corr_help_mentions_ora_only_default():
    parser = build_parser()
    corr_parser = parser._subparsers._group_actions[0].choices["corr"]
    help_text = corr_parser.format_help()
    assert "optionally run GSEA and/or ORA" in help_text
    assert "ora-p-threshold" in help_text
    assert "Force GSEA on." in help_text
    assert "want both analyses." in help_text
    assert "--no-gsea" in help_text
