import gseapy
import pandas as pd


def gsea(gene_names, gene_scores, gene_set=None, out_dir=None, fmt="png",
         no_plot=False):
    """Perform pre-rank GSEA on a set of gene labels and scores."""
    # Create a dataframe with the gene names and scores

    rnk = pd.DataFrame(list(zip(gene_names, gene_scores)))
    gene_set = gene_set
    # Perform GSEA
    gsea_results = gseapy.prerank(rnk=rnk,
                                  gene_set=gene_set,
                                  outdir=out_dir,
                                  permutation_num=0,
                                  format=fmt,
                                  no_plot=no_plot)
    return gsea_results
