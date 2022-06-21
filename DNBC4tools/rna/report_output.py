import os,argparse,shutil,glob
import scipy.io
from scipy.sparse import csr_matrix
import anndata
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i','--indir', metavar='FILE',default=os.getcwd())
parser.add_argument('--no_bam', action='store_true')
args = parser.parse_args()
indir = args.indir

output = os.path.join(indir,'output')
os.makedirs('%s/output/attachment'%indir, exist_ok=True)
def summary_file(indir,output):
    if os.path.exists('%s/04.report/metrics_summary.xls'%indir):
        shutil.copy("%s/04.report/metrics_summary.xls"%indir,'%s/output'%indir)
    if os.path.exists("%s/01.data/raw_matrix"%indir):
        if os.path.exists('%s/raw_matrix'%output):
            shutil.rmtree('%s/raw_matrix'%output)
        shutil.copytree("%s/01.data/raw_matrix"%indir,'%s/raw_matrix'%output,dirs_exist_ok=True)
    if os.path.exists("%s/02.count/filter_matrix"%indir):
        if os.path.exists('%s/filter_matrix'%output):
            shutil.rmtree('%s/filter_matrix'%output)
        shutil.copytree("%s/02.count/filter_matrix"%indir,'%s/filter_matrix'%output,dirs_exist_ok=True)    
    if os.path.exists('%s/03.analysis/Clustering/clustering_plot.png'%indir):
        if os.path.exists('%s/attachment/Clustering'%output):
            shutil.rmtree('%s/attachment/Clustering'%output)
        shutil.copytree("%s/03.analysis/Clustering"%indir,'%s/attachment/Clustering'%output,dirs_exist_ok=True)
        os.remove("%s/attachment/Clustering/cell_report.csv"%output)
        os.remove("%s/attachment/Clustering/cluster_cell.stat"%output)
    if os.path.exists('%s/03.analysis/QC/raw_QCplot.png'%indir):
        if os.path.exists('%s/attachment/QC'%output):
            shutil.rmtree('%s/attachment/QC'%output)
        shutil.copytree("%s/03.analysis/QC"%indir,'%s/attachment/QC'%output,dirs_exist_ok=True)

def ReadPISA(path):
    mat = scipy.io.mmread(path+"/"+"matrix.mtx.gz").astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat,dtype="float32")
    genes = pd.read_csv(path+'/'+'features.tsv.gz', header=None, sep='\t')
    var_names = genes[0].values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    adata.obs_names = pd.read_csv(path+'/'+'barcodes.tsv.gz', header=None)[0].values
    adata.var_names_make_unique()
    return adata

def glob_file(indir):
    html_list = glob.glob('%s/04.report/*.html'%indir)
    if html_list:
        shutil.copy("%s"%html_list[0],'%s/output'%indir)

def main():
    if args.no_bam:
        pass
    else:
        if os.path.exists('%s/02.count/anno_decon_sorted.bam'%indir):
            if(os.path.exists("%s/anno_decon_sorted.bam"%output)):
                os.remove("%s/anno_decon_sorted.bam"%output)
            shutil.move("%s/02.count/anno_decon_sorted.bam"%indir,'%s'%output)
        if os.path.exists('%s/02.count/anno_decon_sorted.bam.bai'%indir):
            if(os.path.exists("%s/anno_decon_sorted.bam.bai"%output)):
                os.remove("%s/anno_decon_sorted.bam.bai"%output)
            shutil.move("%s/02.count/anno_decon_sorted.bam.bai"%indir,'%s'%output)
    summary_file(indir,output)
    adata = ReadPISA('%s/02.count/filter_matrix'%indir)
    adata.write("%s/output/filter_feature.h5ad"%indir)
    glob_file(indir)

if __name__ == '__main__':
    main()

    
