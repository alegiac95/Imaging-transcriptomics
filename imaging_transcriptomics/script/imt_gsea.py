#!/usr/bin/env python3

import pickle
import argparse
from pathlib import Path
import gseapy


def parse_args():
    DESCRIPTION = ""
    EPILOG = ""
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG)
    parser.add_argument("-i", "--input", type=str, required=False,
                        help="Path to the input. This MUST be a pickle file, "
                             "(i.e. a file ending in .pkl), created by "
                             "running an imaging transcriptomics "
                             "analysis with the imagingtranscriptomics "
                             "script.")
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="Path where the results will be saved. If this "
                             "is not provided the same path as the input "
                             "will be used.")
    parser.add_argument("-g", "--geneset", type=str, required=False,
                        default="lake",
                        help="Name of the gene set to use. Some of the "
                             "avilable gene sets are:\n"
                             "- 'lake': \n"
                             "- 'pooled' \n"
                             "- 'kegg' \n"
                             "The 'lake' and 'pooled' gene sets are inluded "
                             "in the imaging transcriptomics package while "
                             "all the other gene sets are available in the "
                             "gseapy package. If you want to see all the "
                             "available gene sets, please run the "
                             "'imt_gsea -g avail' command.")
    parser.add_argument("-m", "--max_genes", type=int, required=False,
                        default=500,
                        help="Maximum number of genes to use in the "
                             "analysis. Default is 500.")

    return parser.parse_args()


def main():
    parsed = parse_args()
    if parsed.geneset == "avail":
        avail_gene_sets = ["lake", "pooled"] + gseapy.get_library_name()
        print("The following gene sets are available:")
        for gene_set in avail_gene_sets:
            print(f"- {gene_set}")
    else:
        geneset = parsed.geneset
        if parsed.input is None:
            raise ValueError("Please provide an input file.")
        infile = Path(parsed.input)
        if not infile.suffix == ".pkl":
            raise ValueError("The input file must be a pickle file.")
        outdir = Path(parsed.output) if parsed.output else infile.parent
        with open(infile, "rb") as f:
            transcriptomics = pickle.load(f)
        transcriptomics.gsea(outdir=outdir, gene_set=geneset,
                             gene_limit=parsed.max_genes)


if __name__ == '__main__':
    main()

