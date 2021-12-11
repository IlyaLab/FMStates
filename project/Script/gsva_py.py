"""Execute bioconductors GSVA transformation of gene expression into pathway enrichment.

This python package is modified from GSVA python library (https://pypi.org/project/GSVA/).

Find the official R package here:

https://doi.org/doi:10.18129/B9.bioc.GSVA

And if you find this useful, cite the authors publication:

Hänzelmann S, Castelo R and Guinney J (2013). "GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, pp. 7. doi: 10.1186/1471-2105-14-7, http://www.biomedcentral.com/1471-2105/14/7.

"""
import sys, os
import pandas as pd 
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def gsva_py(expression_df,geneset_df=None,
         method='gsva',
         kcdf='Gaussian',
         abs_ranking=False,
         min_sz=1,
         max_sz=None,
         parallel_sz=0,
         #parallel_type="SOCK",
         mx_diff=True,
         tau=None,
         ssgsea_norm=True,
         verbose=False,
         tempdir= None
         ):

    """GSVA function for use with pandas DataFrame objects

    :param expression_df: REQUIRED: Expression data indexed on gene names column labels as sample ids
    :type expression_df: pandas.DataFrame
    :param geneset_df: REQUIRED: Genesets and their members in a dataframe
    :type geneset_df: pandas.DataFrame
    :param method: Method to employ in the estimation of gene-set enrichment scores per sample. By default this is set to gsva (Hänzelmann et al, 2013) and other options 6 gsva are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores over the samples and, in the case of zscore, it combines them together as their sum divided by the square-root of the size of the gene set, while in the case of plage they are used to calculate the singular value decomposition (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector as pathway activity profile.
    :type method: string Default: 'gsva'   
    :param kcdf: Character string denoting the kernel to use during the non-parametric estimation of the cumulative distribution function of expression levels across samples when method="gsva". By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson". This argument supersedes arguments rnaseq and kernel, which are deprecated and will be removed in the next release.
    :type kcdf: string Default: 'Gaussian'
    :param abs_ranking: Flag used only when mx_diff=TRUE. When abs_ranking=FALSE [default] a modified Kuiper statistic is used to calculate enrichment scores, taking the magnitude difference between the largest positive and negative random walk deviations. When abs.ranking=TRUE the original Kuiper statistic that sums the largest positive and negative random walk deviations, is used. In this latter case, gene sets with genes enriched on either extreme (high or low) will be regarded as 'highly’ activated.
    :type abs_ranking: bool Default: False
    :param min_sz: Minimum size of the resulting gene sets.
    :type min_sz: int Default: 1
    :param max_sz: Maximum size of the resulting gene sets. Leave unset for no limit.
    :type max_sz: int Default: Inf
    :param parallel_sz: Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel_sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.
    :type parallel_sz: int Default: 0
    :param parallel_type: Type of cluster architecture when using snow.
    :type parallel_type: string Default: "SOCK"   
    :param mx_diff: Offers two approaches to calculate the enrichment statistic (ES) from the KS random walk statistic. mx_diff=FALSE: ES is calculated as the maximum distance of the random walk from 0. mx_diff=TRUE (default): ES is calculated as the magnitude difference between the largest positive and negative random walk deviations.
    :type mx_diff: bool Default: True    
    :param tau: Exponent defining the weight of the tail in the random walk performed by both the gsva (Hänzelmann et al., 2013) and the ssgsea (Barbie et al., 2009) methods. By default, this tau=1 when method="gsva" and tau=0.25 when method="ssgsea" just as specified by Barbie et al. (2009) where this parameter is called alpha. Leave unset for defaults.
    :type tau: float    
    :param ssgsea_norm: Logical, set to TRUE (default) with method="ssgsea" runs the SSGSEA method from Barbie et al. (2009) normalizing the scores by the absolute difference between the minimum and the maximum, as described in their paper. When ssgsea_norm=FALSE this last normalization step is skipped.
    :type ssgsea_norm: bool Default: True    
    :param verbose: Gives information about each calculation step.
    :type verbose: bool Default: False
    :param tempdir: Location to write temporary files
    :type tempdir: string Default: System Default
    :returns: pandas.DataFrame
    """
    df = expression_df
    gmt_df = geneset_df

    if not tempdir:
        tempdir =  mkdtemp(prefix="weirathe.",dir=gettempdir().rstrip('/'))
    if verbose:
        sys.stderr.write("Caching to "+tempdir+"\n")
    # Remove genes from the genesets that do not occur in the dataset
    members = gmt_df['member'].unique()
    missing = set(members)-set(df.index)
    original = df.index
    if len(missing) > 0:
        if verbose: sys.stderr.write("WARNING removing "+str(len(missing))+\
          " genes from gene sets that don't exist in the data\n"+\
          ",".join(sorted(list(missing)))+"\n")
    gmt_df = gmt_df[~gmt_df['member'].isin(list(missing))]
    # Write our gene sets
    gmt_df = gmt_df.groupby(['name']).\
        apply(lambda x: "\t".join(sorted(list(x['member'])))).reset_index().rename(columns={0:'members'})
    of = open(os.path.join(tempdir,"gs.gmt"),'w')
    for row in gmt_df.itertuples():
        name = row.name
        description = 'description'
        fields = row.members
        of.write(name+"\t"+description+"\t"+fields+"\n")
    of.close()
    df.to_csv(os.path.join(tempdir,"expr.csv"))
    cur = os.path.dirname(os.path.realpath(__file__))
    rscript = os.path.join(cur,"gsva.r")
    cmd = ["Rscript",rscript]+[str(x) for x in [
    method,kcdf,abs_ranking,min_sz,max_sz,parallel_sz,
    mx_diff,ssgsea_norm,verbose,tempdir]]
    if verbose: sys.stderr.write(" ".join(cmd)+"\n")
    destination = PIPE
    if verbose: destination = sys.stderr
    sp = Popen(cmd,stdout=PIPE,stderr=destination)
    sp.communicate()
    output = pd.read_csv(os.path.join(tempdir,"pathways.csv"),index_col=0)
    output.index = output.index.astype(str)
    output.columns = output.columns.astype(str)
    output.index.name = 'name'
    return output

