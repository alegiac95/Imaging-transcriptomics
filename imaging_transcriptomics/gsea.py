import gseapy
import pandas as pd
from pathlib import Path


def gsea(gene_names, gene_scores, gene_set: str = None, out_dir=None,
         fmt="png",
         no_plot=False):
    """Perform pre-rank GSEA on a set of gene labels and scores."""
    # Create a dataframe with the gene names and scores

    rnk = pd.DataFrame(list(zip(gene_names, gene_scores)))
    if gene_set.lower() == "lake":
        gene_set = Path(__file__).parent.parent / "data" / "geneset_LAKE.gmt"
    elif gene_set.lower() == "pooled":
        gene_set = Path(__file__).parent.parent / "data" / "geneset_Pooled.gmt"
    else:
        gene_set = gene_set
    # Perform GSEA
    gsea_results = gseapy.prerank(rnk=rnk,
                                  gene_set=gene_set,
                                  outdir=out_dir,
                                  permutation_num=0,
                                  format=fmt,
                                  no_plot=no_plot)
    return gsea_results


def ensemble_gsea(gene_names, gene_scores, permuted_gene_scores):
    """ Perform GSEA on all permutations of the genes and get the pvalue
    from the results.
    """


