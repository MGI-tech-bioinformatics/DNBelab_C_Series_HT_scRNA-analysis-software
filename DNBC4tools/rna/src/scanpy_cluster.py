import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
from typing import Tuple, List
import pandas as pd
import numpy as np
import seaborn as sns
import math
from anndata import AnnData
from sklearn.metrics import pairwise_distances
from dnbc4tools.__init__ import __root_dir__

class CellDataAnalyzer:
    """
    A class for single-cell analysis workflows.
    
    Functions:
        _pp_qc: quality control (QC) filtering
        _pp_basicfilter: basic filtering based on QC metrics
        _pp_scanpy_doubletdetect: doublet detection using scrublet
        pp_hvgs
        pp_reduce
        pp_cluster
        pp_deg
        pp_autoanno
    """

    def __init__(self, adata: sc.AnnData):
        """
        Initializes the CellDataAnalyzer object.

        Parameters:
            adata (AnnData): An annotated data matrix of shape (n_obs, n_vars).
        """
        self.adata = adata
    
    
    def _pp_qc(self, mtgene: list = None):
        """
        Performs quality control (QC) filtering based on mitochondrial genes.

        Parameters:
            mtgene (list): A list of mitochondrial gene names.
        """
        is_mt = self.adata.var_names.str.startswith('MT-', 'mt-')
        if mtgene:
            is_mt = is_mt | self.adata.var_names.isin(mtgene)
        self.adata.var['mt'] = is_mt  
        sc.pp.calculate_qc_metrics(
            self.adata, 
            qc_vars=['mt'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )

    def _pp_basicfilter(self,
                        min_cells: int = 3,
                        pct_counts_mt: bool = True,
                        filter_maxgene: bool = True,
                        filter_mingene: bool = True
                        ):
        """
        Performs basic filtering of cells and genes based on QC metrics.

        Parameters:
            min_cells (int): The minimum number of cells that a gene must be detected in.
            pct_counts_mt (bool): Whether to filter cells based on the percentage of mitochondrial genes.
            filter_maxgene (bool): Whether to filter cells based on the maximum number of detected genes.
            filter_mingene (bool): Whether to filter cells based on the minimum number of detected genes.
        """
        if pct_counts_mt:
            threshold = math.ceil(self.adata.obs.pct_counts_mt.quantile(0.95))
            threshold = threshold if threshold >= 5 else 5
            self.adata = self.adata[self.adata.obs.pct_counts_mt <= threshold, :]
        if filter_maxgene:
            threshold = math.ceil(self.adata.obs.n_genes_by_counts.quantile(0.95))
            self.adata = self.adata[self.adata.obs.n_genes_by_counts <= threshold, :]
        if filter_mingene:
            threshold = math.ceil(self.adata.obs.n_genes_by_counts.quantile(0.05))
            threshold = threshold if threshold <= 200 else 200
            self.adata = self.adata[self.adata.obs.n_genes_by_counts >= threshold, :]
        if min_cells:
            sc.pp.filter_genes(self.adata, min_cells=min_cells)
        
    def _pp_scanpy_doubletdetect(self,
                                doublet_threshold: float = None,
                                doublet_rate: float = 0.05
                                ):
        """
        Detect doublets using the scrublet algorithm from scanpy.external.pp module.

        Args:
            doublet_threshold (float, optional): The scrublet score threshold used to classify cells as doublets. Defaults to None.
            doublet_rate (float, optional): Expected doublet rate. Defaults to 0.05.
        """
        scrub_doublet = sc.external.pp.scrublet_simulate_doublets(
            self.adata,
            sim_doublet_ratio=2.0,
            synthetic_doublet_umi_subsampling=1.0,
            random_seed=123
        )
        try:
            sc.external.pp.scrublet(
                self.adata,
                scrub_doublet,
                expected_doublet_rate=doublet_rate,
                threshold=doublet_threshold,
                copy=False,
                verbose=True,
                random_state=123
            )
        except Exception as e:
            cell_num = len(self.adata.obs.index)
            n_comps = cell_num - 1
            sc.external.pp.scrublet(
                self.adata,
                scrub_doublet,
                expected_doublet_rate=doublet_rate,
                threshold=doublet_threshold,
                copy=False,
                verbose=True,
                random_state=123,
                n_prin_comps=n_comps
                )
        self.adata = self.adata[self.adata.obs['predicted_doublet'] != True]


    def _pp_hvgs(self,
                    min_mean: float = 0.0125,
                    max_mean: float = 3,
                    min_disp: float = 0.5,
                    hvgs: int = 2000,
                    normalization: str = 'LogNormalize',
                    flavor: str = 'seurat_v3'
                ):
        """
        Preprocesses data to extract highly variable genes for scRNA-seq analysis.

        min_mean: Minimum mean expression value for a gene to be considered for analysis.
        max_mean: Maximum mean expression value for a gene to be considered for analysis.
        min_disp: Minimum dispersion value for a gene to be considered for analysis.
        hvgs: Number of highly variable genes to select.
        normalization: Type of normalization to apply to the data.
        flavor: Flavor of highly variable gene selection algorithm to use.
        """
        if normalization == 'LogNormalize':
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
        if flavor == 'seurat_v3':
            try:
                sc.pp.highly_variable_genes(self.adata, flavor='seurat_v3', n_top_genes=hvgs)
            except Exception as e:
                sc.pp.highly_variable_genes(self.adata, flavor='seurat_v3', n_top_genes=hvgs, span=1)
        else:
            sc.pp.highly_variable_genes(self.adata, min_mean=min_mean, max_mean=max_mean,
                                        min_disp=min_disp, flavor=flavor, n_top_genes=hvgs)
        self.adata.raw = self.adata
        sc.pp.scale(self.adata, max_value=10)


    def _pp_reduce(self,
                    n_comps: int = 50,
                    svd_solver: str = 'arpack',
                    use_highly_variable: bool = True
                    ):
        """
        Perform PCA dimensionality reduction on the expression matrix of single-cell data.

        Args:
            n_comps (int): Number of principal components to compute. Default is 50.
            svd_solver (str): SVD solver to use. Must be one of {'arpack', 'randomized'}.
                            Default is 'arpack'.
            use_highly_variable (bool): If True, only use highly variable genes for PCA.
                                        Default is True.
        """
        try:
            sc.tl.pca(self.adata, n_comps=n_comps, svd_solver=svd_solver,
                    use_highly_variable=use_highly_variable, random_state=123)
        except Exception as e:
            cell_num = len(self.adata.obs.index)
            n_comps = cell_num - 1
            sc.tl.pca(self.adata, n_comps=n_comps, svd_solver=svd_solver,
                    use_highly_variable=use_highly_variable, random_state=123)


    def _pp_cluster(self,
                    neighbors: int = 15,
                    n_pcs: int = 30,
                    method: str = 'umap',
                    spread: float = 1.0,
                    min_dist: float = 0.25,
                    perplexity: int = 30,
                    learning_rate: int = 1000,
                    cluster_method: str = 'louvain',
                    resolution: int = 0.8
                    ):
        """
        Perform clustering analysis on preprocessed single-cell RNA sequencing data.

        Args:
            neighbors (int): Number of neighbors for construction of the nearest-neighbors graph.
            n_pcs (int): Number of principal components to use.
            method (str): Method used for dimension reduction. Default is "umap".
            spread (float): The effective scale of embedded points. Default is 1.0.
            min_dist (float): The effective minimum distance between embedded points. Default is 0.25.
            perplexity (int): The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. 
                            Larger datasets usually require a larger perplexity. Default is 30.
            learning_rate (int): The learning rate is a parameter that controls the step size at each iteration of gradient descent. Default is 1000.
            cluster_method (str): Clustering method to use. Default is "louvain".
            resolution (int): A parameter value controlling the coarseness of the clustering. Higher values result in more clusters. Default is 0.8.
        """
        try:
            sc.pp.neighbors(self.adata, n_neighbors=neighbors, n_pcs=n_pcs, method='umap', metric='euclidean', random_state=123)
        except Exception as e:
            n_comp = len(self.adata.uns['pca']['variance'])
            sc.pp.neighbors(self.adata, n_neighbors=neighbors, n_pcs=n_comp, method='umap', metric='euclidean', random_state=123)

        if method == 'umap':
            sc.tl.umap(self.adata, min_dist=min_dist, spread=spread, random_state=123)
        if method == 'tsne':
            sc.tl.tsne(self.adata, perplexity=perplexity, learning_rate=learning_rate, random_state=123)

        if cluster_method == 'leiden':
            sc.tl.leiden(self.adata, resolution=resolution, random_state=123)
        if cluster_method == 'louvain':
            sc.tl.louvain(self.adata, resolution=resolution, random_state=123)


    def _pp_deg(self,
                groupby: str = 'leiden/louvain',
                method: str = 'wilcoxon',
                n_genes: int = 200,
                ):
        """
        Differential expression analysis using `sc.tl.rank_genes_groups` function from `scanpy` package.
        
        Args:
            groupby: The categorical annotation grouping to consider. If the dataset is clustered
                (i.e. clusters are defined), then `groupby` defaults to `'leiden'`, otherwise `'louvain'`.
            method: The method used to calculate gene-level summary statistics. Default is `'wilcoxon'`.
            n_genes: The number of marker genes to keep. Default is `200`.
        """
  
        self.adata.uns['log1p'] = {'base': None}
        groupby = 'louvain' if 'louvain' in self.adata.obs_keys() else 'leiden'

        if len(self.adata.obs[groupby].cat.categories) > 1 :
            sc.tl.rank_genes_groups(self.adata, groupby, pts = True,
                                method=method, n_genes=n_genes)
        else:
            pass

    # def _pp_autoanno(self,
    #                 species: str = None
    #                 ):
    #     celltypist.models.models_path = '%s/config/cellAnno'%__root_dir__
    #     if species == 'Human':
    #         modeltype = 'Immune_All_Low.pkl'
    #     elif species == 'Mouse':
    #         modeltype = 'Adult_Mouse_Gut.pkl'
    #     else:
    #         modeltype = None
        
    #     predictions = celltypist.annotate(
    #         self.adata, model = modeltype , majority_voting=True, over_clustering = self.adata.obs['louvain'].astype(str))
    #     self.adata= predictions.to_adata()

class AnnoMCA_HCL():

    def __init__(self, reference_df: pd.DataFrame, query_adata: AnnData, n_cores: int):
        """
        Initialize the AnnoMCA_HCL class.

        Args:
            reference_df (pd.DataFrame): A pandas dataframe containing reference expression data.
            query_adata (anndata.AnnData): An Anndata object containing query expression data.
            n_cores (int): The number of cores to use for processing. 
        """
        self.reference_df = reference_df
        self.query_adata = query_adata
        self.n_cores = n_cores

    def adata_to_query(self) -> pd.DataFrame:
        """
        Convert the query expression data in Anndata format to a pandas dataframe.

        Returns:
            pd.DataFrame: A pandas dataframe containing the query expression data.
        """
        query_df = pd.DataFrame(self.query_adata.raw.X.A,
                                index=self.query_adata.obs.index,
                                columns=self.query_adata.raw.var.index)
        return query_df

    def process_query_df(self, query_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process the query dataframe by subseting and log-transforming the data.

        Args:
            query_df (pd.DataFrame): A pandas dataframe containing the query expression data.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two pandas dataframes: 
                1. A subset of the query dataframe containing only the genes present in the reference dataset, 
                2. A log-transformed version of the reference dataframe.
        """
        tst_expr_subset = []
        unfound_genes = []
        found_genes = []
        for gene in self.reference_df.columns:
            if gene in query_df.columns:
                found_genes.append(gene)
                tst_expr_subset.append(query_df[gene])

            else:
                unfound_genes.append(gene)
                _t = pd.Series(0, index=query_df.index, name=gene)
                tst_expr_subset.append(_t)

        tst_expr_subset = pd.DataFrame(tst_expr_subset).T
        tst_expr_subset = 1e5 * tst_expr_subset / tst_expr_subset.values.sum(1, keepdims=True)
        tst_expr_subset = np.log(tst_expr_subset + 1)
        ref_expr_trans = np.log(self.reference_df + 1)
        return tst_expr_subset, ref_expr_trans
    
    def pearson_matrix(self, tst_expr_subset: pd.DataFrame, ref_expr_trans: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate the Pearson correlation matrix between the query and reference data.

        Args:
            tst_expr_subset (pd.DataFrame): A subset of the query dataframe containing only the genes present in the reference dataset.
            ref_expr_trans (pd.DataFrame): A log-transformed version of the reference dataframe.

        Returns:
            pd.DataFrame: A pandas dataframe containing the Pearson correlation matrix.
        """
        assert np.all(tst_expr_subset.columns == ref_expr_trans.columns)
        pear_matrix = 1 - pairwise_distances(tst_expr_subset, ref_expr_trans, metric='correlation', n_jobs=self.n_cores)
        pear_matrix = pd.DataFrame(pear_matrix, index=tst_expr_subset.index, columns=ref_expr_trans.index)

        closest_ref = np.argmax(pear_matrix.values, axis=1)
        closest_ref_name = pear_matrix.columns.values[closest_ref]
        score = pear_matrix.values[np.arange(len(closest_ref)), closest_ref]
        celltype = pd.Series(closest_ref_name, index=pear_matrix.index, name='celltype')
        celltype_probility = pd.Series(score, index=pear_matrix.index, name='celltype_probility')
        celltype_df = pd.DataFrame({'celltype_probility': celltype_probility, 'celltype': celltype})
        return celltype_df
    
def run_AnnoMCA_HCL(reference_df: pd.DataFrame, query_adata: AnnData, n_cores: int) -> AnnData:
    """
    Run AnnoMCA_HCL algorithm to annotate query cells based on reference data.

    Args:
    - reference_df: A pandas DataFrame containing the expression data of the reference cells.
    - query_adata: An AnnData object containing the expression data of the query cells.
    - n_cores: An integer specifying the number of cores to use for parallel computing.

    Returns:
    - An AnnData object with the annotated cell types added as a new column to the .obs dataframe.
    """
    anno_model = AnnoMCA_HCL(reference_df, query_adata, n_cores)
    query_df = anno_model.adata_to_query()
    tst_expr_subset, ref_expr_trans = anno_model.process_query_df(query_df)
    celltype_df = anno_model.pearson_matrix(tst_expr_subset, ref_expr_trans)
    merged_df = pd.concat([query_adata.obs, celltype_df], axis=1)
    for category, group in merged_df.groupby('louvain'):
        most_common_value = group['celltype'].value_counts().idxmax()
        merged_df.loc[group.index, 'celltype'] = most_common_value
    query_adata.obs = merged_df
    return query_adata


def mtgene_file_read(mtgenefile: str) -> List[str]:
    """
    Read a file containing mitochondrial gene names and return a list of gene names.

    Args:
    - mtgenefile: A string specifying the path of the file containing mitochondrial gene names.

    Returns:
    - A list of mitochondrial gene names.
    """
    mtgene = []
    df = pd.read_csv(mtgenefile, header=None)
    for line in list(df.iloc[:, 0]):
        line = line.strip()
        mtgene.append(line)
    return mtgene


def get_markers(adata: AnnData, cutoff: float = 0.05) -> pd.DataFrame:
    """
    Obtain differentially expressed genes that define each cluster in the AnnData object.
    
    Parameters:
        adata (AnnData): AnnData object containing gene expression data and cluster labels.
        cutoff (float, optional): Adjusted p-value cutoff for differential expression. Default is 0.05.
    
    Returns:
        pd.DataFrame: Table containing differentially expressed genes for each cluster, along with their scores,
                      average log2 fold change, p-values, and percentages in each group.
    """
    markertable = sc.get.rank_genes_groups_df(adata,group=None)
    markertable = markertable[~pd.isnull(markertable['logfoldchanges'])]
    markertable.columns = ['cluster','gene','score','avg_log2FC','p_val','p_val_adj','pct.1','pct.2']
    markertable = markertable.round({'avg_log2FC':3,'pct.1':3,'pct.2':3})
    markertable = markertable.loc[markertable['p_val_adj'] < cutoff, ]
    markertable = markertable.loc[markertable['avg_log2FC'] > 0].sort_values('cluster', ascending=True).groupby('cluster').head(30)
    return markertable

def get_cluster(adata: AnnData, cluster_method: str = 'louvain') -> pd.DataFrame:
    """
    Obtain cluster information from the AnnData object and assign unique names to each cluster based on the cluster
    method used for clustering and the number of cells in each cluster. If cell type annotation information is available,
    predict the cell types for each cluster based on the majority voting of cell types for cells in each cluster.
    
    Parameters:
        adata (AnnData): AnnData object containing gene expression data and cluster labels.
        cluster_method (str, optional): Method used for clustering. Default is 'louvain'.
    
    Returns:
        pd.DataFrame: Table containing cluster information, including number of genes and UMIs per cell, cluster IDs,
                      UMAP coordinates, unique cluster names, predicted number of cells for each cell type (if available),
                      and predicted cell type (if available).
    """
    umap_matrix = adata.obsm.to_df()[['X_umap1','X_umap2']]
    # if 'majority_voting' in adata.obs:
    #     obs_data = adata.obs[['n_genes_by_counts','total_counts','%s'%cluster_method,'majority_voting']]
    if 'celltype' in adata.obs:
        obs_data = adata.obs[['n_genes_by_counts','total_counts','%s'%cluster_method,'celltype']]
    else:
        obs_data = adata.obs[['n_genes_by_counts','total_counts','%s'%cluster_method]]
    
    merge_data = pd.concat([obs_data,umap_matrix],axis=1)
    merge_data['Cluster_number'] =merge_data['%s'%cluster_method].map(merge_data['%s'%cluster_method].value_counts(ascending=False))
    merge_data['Cluster'] = merge_data['%s'%cluster_method].astype(str) + ' CellsNum: ' + merge_data['Cluster_number'].astype(str)
    merge_data = merge_data.drop(['Cluster_number'],axis=1)
    # if 'majority_voting' in merge_data:
    #     merge_data['Predict_number'] =merge_data['majority_voting'].map(merge_data['majority_voting'].value_counts(ascending=False))
    if 'celltype' in merge_data:
        merge_data['Predict_number'] =merge_data['celltype'].map(merge_data['celltype'].value_counts(ascending=False))
        merge_data['Predicted cell type'] = merge_data['celltype'].astype(str) + ' : ' + merge_data['Predict_number'].astype(str)
        merge_data = merge_data.drop(['celltype'],axis=1)
        merge_data.columns = ['nGene','nUMI','%s'%cluster_method,'UMAP_1','UMAP_2','Cluster','Predict_number','Predicted cell type']
    else:
        merge_data.columns = ['nGene','nUMI','%s'%cluster_method,'UMAP_1','UMAP_2','Cluster']
    merge_data = merge_data.sort_values(by=['%s'%cluster_method],ascending=True)
    return merge_data


def draw_qcfigure(adata,
                  ):
    params = {'figure.figsize': (5, 5),
            'axes.labelsize': 'medium',
            'figure.dpi': 150,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'xtick.labelsize':'medium',
            'ytick.labelsize':'x-small'}
    plt.rcParams.update(params)
    meta = adata.obs

    threshold1 = math.ceil(meta['pct_counts_mt'].quantile(0.995))
    threshold2 = math.ceil(meta['n_genes_by_counts'].quantile(0.995))
    threshold3 = math.ceil(meta['total_counts'].quantile(0.995))
    meta_threshold1 = meta[meta['pct_counts_mt'] <= threshold1]
    meta_threshold2 = meta[meta['n_genes_by_counts'] <= threshold2]
    meta_threshold3 = meta[meta['total_counts'] <= threshold3]

    boxprops = dict(linewidth=1, edgecolor='black',zorder=1)
    whiskerprops = dict(linewidth=1, color='black')
    medianprops = dict(linewidth=1, color='black')
    if meta['pct_counts_mt'].sum() > 0:
        fig, ax = plt.subplots(1,3)
        sns.violinplot(
            y="pct_counts_mt", 
            data=meta_threshold1, 
            color="#7570B3",
            inner=None,
            bw=.2, 
            saturation=1,
            scale="width",
            ax=ax[2],
            cut=0,
            linewidth=0.5
            )
        sns.boxplot(
            y="pct_counts_mt", 
            data=meta, 
            ax=ax[2], 
            width=0.2, 
            color="white",
            fliersize=0,
            showfliers=False,
            boxprops=boxprops, 
            whiskerprops=whiskerprops,
            medianprops=medianprops,
            showcaps=False)
        ax[2].axes.set_ylabel('')
        ax[2].axes.set_xlabel('mito.percent')
    else:
        fig, ax = plt.subplots(1,2)
    sns.violinplot(
        y="n_genes_by_counts", 
        data=meta_threshold2 , 
        color="#1B9E77",
        inner=None,
        bw=.2, 
        saturation=1,
        scale="width",
        ax=ax[0],
        cut=0,
        linewidth=0.5
        )
    sns.boxplot(
        y="n_genes_by_counts", 
        data=meta, 
        ax=ax[0], 
        width=0.2, 
        color="white",
        fliersize=0,
        showfliers=False,
        boxprops=boxprops, 
        whiskerprops=whiskerprops,
        medianprops=medianprops,
        showcaps=False
        )
    ax[0].axes.set_ylabel('')
    ax[0].axes.set_xlabel('genes')

    sns.violinplot(
        y="total_counts", 
        data=meta_threshold3, 
        color="#D95F02",
        inner=None,
        bw=.2, 
        saturation=1,
        scale="width",
        ax=ax[1],
        cut=0,
        linewidth=0.5
        )
    sns.boxplot(
        y="total_counts", 
        data=meta,
        ax=ax[1], 
        width=0.2, 
        color="white",
        fliersize=0,
        showfliers=False,
        boxprops=boxprops, 
        whiskerprops=whiskerprops,
        medianprops=medianprops,
        showcaps=False)
    ax[1].axes.set_ylabel('')
    ax[1].axes.set_xlabel('counts')

    fig.tight_layout()
    return fig
