#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from imaging_transcriptomics import atlas_table, run_corr, run_pls


def _default_output_dir(input_path: str, method: str) -> Path:
    path = Path(input_path)
    name = path.name
    if name.endswith(".nii.gz"):
        stem = name[:-7]
    else:
        stem = path.stem
    return path.parent / f"Imt_{stem}_{method}"


def parse_cmdline():
    description = "Run imaging transcriptomics 2.0 analyses with atlas-aware scan extraction."
    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(dest="command", required=True)

    atlas_parser = subparsers.add_parser("atlases", help="List packaged and buildable atlas presets.")
    atlas_parser.add_argument(
        "--packaged-only",
        action="store_true",
        help="Show only atlases that ship ready-to-run expression data in this branch.",
    )

    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument("-i", "--input", required=True, help="Input scan, vector file, or volumetric image.")
    shared.add_argument(
        "--input-rh",
        default=None,
        help="Right-hemisphere surface file when using surface inputs.",
    )
    shared.add_argument("-o", "--output", default=None, help="Output directory. Defaults to a method-specific folder next to the input.")
    shared.add_argument(
        "-a",
        "--atlas",
        default="dk",
        help="Atlas preset to use. Examples: dk, schaefer-100, schaefer-200, schaefer-400, destrieux, glasser-360.",
    )
    shared.add_argument(
        "--hemisphere",
        choices=["left", "both"],
        default="left",
        help="Use left-hemisphere only or both hemispheres from the packaged abagen matrix.",
    )
    shared.add_argument(
        "-r",
        "--regions",
        choices=["all", "cort+sub", "cort"],
        default="all",
        help="Region scope to analyze.",
    )
    shared.add_argument(
        "--space",
        default=None,
        help="Source space for the input data, e.g. MNI152, fsaverage, fsLR, CIVET.",
    )
    shared.add_argument(
        "-p",
        "--permutations",
        type=int,
        default=1000,
        help="Number of permutations or spins to use.",
    )
    shared.add_argument(
        "--seed",
        type=int,
        default=1234,
        help="Random seed for permutation and null-model generation.",
    )
    shared.add_argument(
        "--null-method",
        choices=["auto", "vasa", "alexander_bloch", "moran", "random"],
        default="auto",
        help="Spatial null model for cortical permutations. 'auto' prefers Vasa spins and falls back to random within-hemisphere shuffles.",
    )
    shared.add_argument(
        "--geneset",
        default="lake",
        help="Gene set library to use when GSEA is enabled.",
    )
    shared.add_argument(
        "--no-gsea",
        action="store_true",
        help="Skip GSEA and only write the regional and gene tables plus plots.",
    )

    corr_parser = subparsers.add_parser("corr", parents=[shared], help="Run correlation-based imaging transcriptomics.")
    corr_parser.set_defaults(method="corr")

    pls_parser = subparsers.add_parser("pls", parents=[shared], help="Run PLS-based imaging transcriptomics.")
    pls_group = pls_parser.add_mutually_exclusive_group(required=True)
    pls_group.add_argument("--ncomp", type=int, default=None, help="Number of PLS components to retain.")
    pls_group.add_argument("--var", type=float, default=None, help="Cumulative variance target for selecting PLS components.")
    pls_parser.set_defaults(method="pls")
    return parser.parse_args()


def main():
    parsed = parse_cmdline()
    if parsed.command == "atlases":
        print(atlas_table(packaged_only=parsed.packaged_only).to_string(index=False))
        return

    output_dir = Path(parsed.output) if parsed.output else _default_output_dir(parsed.input, parsed.method)
    run_gsea = not parsed.no_gsea

    if parsed.method == "corr":
        run_corr(
            parsed.input,
            atlas=parsed.atlas,
            hemisphere=parsed.hemisphere,
            regions=parsed.regions,
            source_space=parsed.space,
            input_rh=parsed.input_rh,
            n_permutations=parsed.permutations,
            null_method=parsed.null_method,
            output_dir=output_dir,
            run_gsea=run_gsea,
            gene_set=parsed.geneset,
            seed=parsed.seed,
        )
        return

    run_pls(
        parsed.input,
        atlas=parsed.atlas,
        hemisphere=parsed.hemisphere,
        regions=parsed.regions,
        source_space=parsed.space,
        input_rh=parsed.input_rh,
        n_components=parsed.ncomp,
        var=parsed.var,
        n_permutations=parsed.permutations,
        null_method=parsed.null_method,
        output_dir=output_dir,
        run_gsea=run_gsea,
        gene_set=parsed.geneset,
        seed=parsed.seed,
    )


if __name__ == "__main__":
    main()
