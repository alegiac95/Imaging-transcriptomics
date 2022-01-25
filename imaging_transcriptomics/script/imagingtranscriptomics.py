#!/usr/bin/env python3

import argparse
import imaging_transcriptomics as imt
from pathlib import Path


def parse_cmdline():
    DESCRIPTION = """Perform imaging transcriptomics analysis of a 
    neuroimaging data."""
    EPILOG = """
    If you use this software in your work, please cite:
    
    * Imaging transcriptomics: Convergent cellular, transcriptomic, 
    and molecular neuroimaging signatures in the healthy adult human brain.* 
    Daniel Martins, Alessio Giacomel, Steven CR Williams, Federico Turkheimer,
    Ottavia Dipasquale, Mattia Veronese, PET templates working group. Cell 
    Reports; doi: [https://doi.org/10.1016/j.celrep.2021.110173]
    (https://doi.org/10.1016/j.celrep.2021.110173)
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG)
    # IO arguments
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input file, can be a neuroimaging file (e.g., "
                             ".nii[.gz] or a text file (e.g, .txt, .tsv, "
                             ".csv) containing the values in a column")
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="Output directory, if not provided, the same of "
                             "the input file is used")
    # Control arguments
    parser.add_argument("-r", "--regions", type=str, required=False,
                        choices=["all", "cort+sub", "cort"], default="all",
                        help="Regions to be used for the analysis, can be "
                             "either 'all' (default), 'cort+sub' or 'cort'."
                             "The behaviour with 'all' is the same as "
                             "'cort+sub' and will use all regions to perform "
                             "the analysis, while with 'cort' will use only "
                             "the cortical regions.")
    parser.add_argument("--no-gsea", action="store_false", required=False,
                        help="If True perform GSEA analysis, otherwise skip.")
    parser.add_argument("--geneset", type=str, required=False, default="lake",
                        help="geneset to use for the GSEA analysis. Can be "
                             "either 'lake' (default), 'pooled' or any of "
                             "the genesets included in the gseapy package.")
    parser.add_argument("--max_genes", type=int, required=False,
                        default=500,
                        help="Maximum number of genes to use in the "
                             "analysis. Default is 500.")
    subparser = parser.add_subparsers(title="method", dest="method")
    parse_corr = subparser.add_parser("corr")
    parse_corr.add_argument("--cpu", type=int, required=False, default=4,
                            help="Number of CPUs to use for the analysis.")
    parse_pls = subparser.add_parser("pls")
    pls_group = parse_pls.add_mutually_exclusive_group(required=True)
    pls_group.add_argument("--ncomp", type=int, help="Number of "
                                                     "PLS components.")
    pls_group.add_argument("--var", type=float,
                           help="Percentage of variance to extract form "
                                "the data.")
    return parser.parse_args()


def main():
    parsed = parse_cmdline()
    regions = parsed.regions
    gsea = parsed.no_gsea
    geneset = parsed.geneset
    infile = Path(parsed.input)
    input_name = infile.stem
    outdir = Path(parsed.output) if parsed.output else infile.parent
    if parsed.method == "corr":
        if infile.suffix in [".txt", ".tsv", ".csv"]:
            transcriptomics = imt.ImagingTranscriptomics.from_file(
                infile,
                method="corr",
                regions=regions)
        elif str().join(infile.suffixes) in [
            ".nii", ".nii.gz"]:
            transcriptomics = imt.ImagingTranscriptomics.from_scan(
                infile,
                method="corr",
                regions=regions)
        n_cpu = parsed.cpu
    elif parsed.method == "pls":
        pls_arg = {
            "n_components": parsed.ncomp,
            "var": parsed.var
        }
        if infile.suffix in [".txt", ".tsv", ".csv"]:
            transcriptomics = imt.ImagingTranscriptomics.from_file(
                infile,
                method="pls",
                regions=regions, **pls_arg)
        elif str().join(infile.suffixes) in [".nii", ".nii.gz"]:
            transcriptomics = imt.ImagingTranscriptomics.from_scan(
                infile,
                method="pls",
                regions=regions, **pls_arg)
        n_cpu = 4
    else:
        raise ValueError("Method not recognized")
    transcriptomics.run(outdir,
                        scan_name=input_name,
                        gsea=gsea,
                        gene_set=geneset,
                        n_cpu=n_cpu,
                        gene_limit=parsed.max_genes)
    # PLOTTING and PDF creation
    imt.reporting.make_pdf(transcriptomics_data=transcriptomics,
                           save_dir=Path(outdir) / f"Imt_{input_name}_"
                                                   f"{transcriptomics.method}",
                           name=str(input_name),
                           scanname=str(infile.name))


if __name__ == "__main__":
    main()
