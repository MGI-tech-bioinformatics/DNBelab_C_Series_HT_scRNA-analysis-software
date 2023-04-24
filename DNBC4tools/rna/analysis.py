import os
from dnbc4tools.tools.utils import str_mkdir,judgeFilexits,change_path,read_json,rm_temp
from dnbc4tools.__init__ import __root_dir__

class Analysis:
    def __init__(self,args):
        self.name = args.name
        self.genomeDir = args.genomeDir
        self.outdir = os.path.join(args.outdir,args.name)
    
    def run(self):
        judgeFilexits(
            self.genomeDir,
            '%s/02.count/filter_feature.h5ad'%self.outdir,
            )
        str_mkdir('%s/03.analysis'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        genomeDir = os.path.abspath(self.genomeDir)
        indexConfig = read_json('%s/ref.json'%genomeDir)
        species = indexConfig['species']
        mtgenes = indexConfig['mtgenes']

        change_path()
        import matplotlib
        matplotlib.use('Agg')
        import warnings
        import pandas as pd
        warnings.filterwarnings("ignore")
        import matplotlib.pyplot as plt
        import scanpy as sc
        from dnbc4tools.rna.src.scanpy_cluster import CellDataAnalyzer,run_AnnoMCA_HCL,mtgene_file_read,get_markers,get_cluster,draw_qcfigure

        print('\nDimensionality reduction, Clustering.')
        rm_temp('%s/03.analysis/cluster_annotation.png'%self.outdir)

        adata = sc.read_h5ad('%s/02.count/filter_feature.h5ad'%self.outdir)
        scanpy_object = CellDataAnalyzer(adata)
        if mtgenes != 'None':
            mtgenelist = mtgene_file_read(mtgenes)
            scanpy_object._pp_qc(mtgenelist)
        else:
            scanpy_object._pp_qc()
        fig1 = draw_qcfigure(scanpy_object.adata)
        fig1.savefig('%s/03.analysis/raw_QCplot.png'%self.outdir, dpi = 150)
        scanpy_object.adata.obs.to_csv("%s/03.analysis/raw_qc.xls"%self.outdir,index=True,sep='\t')
        scanpy_object._pp_basicfilter()
        
        if len(scanpy_object.adata.obs.index) > 100:
            scanpy_object._pp_scanpy_doubletdetect()
        
        fig2 = draw_qcfigure(scanpy_object.adata)
        fig2.savefig('%s/03.analysis/filter_QCplot.png'%self.outdir, dpi = 150)
        scanpy_object._pp_hvgs()
        scanpy_object._pp_reduce()
        scanpy_object._pp_cluster()
        scanpy_object._pp_deg()

        if 'rank_genes_groups' in scanpy_object.adata.uns:
            marker_table = get_markers(scanpy_object.adata)
            marker_table.to_csv("%s/03.analysis/marker.csv"%self.outdir,index=False)
        else:
            with open("%s/03.analysis/marker.csv"%self.outdir, 'w') as markerfile:
                markerfile.write('cluster,gene,score,avg_log2FC,p_val,p_val_adj,pct.1,pct.2')

        if species in ["human","Human","hg19","hg38","Homo_sapiens"]:
            try:
                # scanpy_object._pp_autoanno("Human")
                # sc.pl.umap(scanpy_object.adata, color = ['majority_voting', 'louvain'], legend_loc = 'on data')
                ref_df = pd.read_csv('%s/config/cellAnno/HCL.ref.csv.gz'%__root_dir__, sep=',', index_col=0).T
                scanpy_object.adata = run_AnnoMCA_HCL(ref_df,scanpy_object.adata,4)
                sc.pl.umap(scanpy_object.adata, color = ['celltype', 'louvain'], legend_loc = 'on data')
                plt.savefig("%s/03.analysis/cluster_annotation.png"%self.outdir)
            except Exception as e:
                sc.pl.umap(scanpy_object.adata, color = ['louvain'], legend_loc = 'on data')
                plt.savefig("%s/03.analysis/cluster.png"%self.outdir)
        elif species in ["mouse","Mouse","mm10","Mus_musculus"]:
            try:
                # scanpy_object._pp_autoanno("Mouse")
                # sc.pl.umap(scanpy_object.adata, color = ['majority_voting', 'louvain'], legend_loc = 'on data')
                ref_df = pd.read_csv('%s/config/cellAnno/MCA.ref.csv.gz'%__root_dir__, sep=',', index_col=0).T
                scanpy_object.adata = run_AnnoMCA_HCL(ref_df,scanpy_object.adata,4)
                sc.pl.umap(scanpy_object.adata, color = ['celltype', 'louvain'], legend_loc = 'on data')
                plt.savefig("%s/03.analysis/cluster_annotation.png"%self.outdir)
            except Exception as e:
                sc.pl.umap(scanpy_object.adata, color = ['louvain'], legend_loc = 'on data')
                plt.savefig("%s/03.analysis/cluster.png"%self.outdir)
        else:
            sc.pl.umap(scanpy_object.adata, color = ['louvain'], legend_loc = 'on data')
            plt.savefig("%s/03.analysis/cluster.png"%self.outdir)
        cluster_table = get_cluster(scanpy_object.adata)
        cluster_table.to_csv("%s/03.analysis/cluster.csv"%self.outdir,index=True)
        scanpy_object.adata.write("%s/03.analysis/QC_Clutser.h5ad"%self.outdir)

def analysis(args):
    Analysis(args).run()

def helpInfo_analysis(parser):
    parser.add_argument(
        '--name',
        metavar='NAME',
        required=True,
        help='Sample name.'
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '--outdir',
        metavar='DIR',
        help='output dir, [default: current directory].',
        default=os.getcwd()
        )
    return parser