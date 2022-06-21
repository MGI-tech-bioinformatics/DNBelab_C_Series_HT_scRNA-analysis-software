import os
from .utils import str_mkdir,logging_call,judgeFilexits,change_path,rscript_path
from DNBC4tools.__init__ import _root_dir

class Analysis:
    def __init__(self,args):
        self.name = args.name
        self.matrix = args.matrix
        self.qcdim = args.qcdim
        self.clusterdim = args.clusterdim
        self.doubletpercentage = args.doubletpercentage
        self.mitpercentage = args.mitpercentage
        self.mtgenes = args.mtgenes
        self.minfeatures = args.minfeatures
        self.PCusage = args.PCusage
        self.resolution = args.resolution
        self.species = args.species
        self.outdir = os.path.join(args.outdir,args.name)
    
    def run(self):
        judgeFilexits(self.matrix)
        str_mkdir('%s/03.analysis/QC'%self.outdir)
        str_mkdir('%s/03.analysis/Clustering'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        new_rscript = rscript_path()

        qc_cmd = '%s %s/rna/QC_analysis.R -I %s -D %s -P %s -M %s -MP %s -F %s -B %s -O %s/03.analysis'\
            %(new_rscript,_root_dir,self.matrix,self.qcdim,self.doubletpercentage,self.mtgenes,self.mitpercentage,self.minfeatures,self.name,self.outdir)
        cluster_cmd = '%s %s/rna/Cluster_analysis.R -I %s/03.analysis/QC -D %s -PC %s -RES %s -O %s/03.analysis -SP %s' \
            %(new_rscript,_root_dir,self.outdir,self.clusterdim,self.PCusage,self.resolution,self.outdir,self.species)
        logging_call(qc_cmd,'analysis',self.outdir)
        logging_call(cluster_cmd,'analysis',self.outdir)
        
def analysis(args):
    Analysis(args).run()

def parse_analysis(parser):
    parser.add_argument('--name',metavar='NAME',required=True,help='Sample name.')
    parser.add_argument('--matrix',metavar='DIR',required=True,help='Count filter matrix dir, contain barcodes.tsv,features.tsv,matrix.mtx.')
    parser.add_argument('--qcdim',metavar='INT',type=int,default=20,help="DoubleFinder's PCs parameter, the number of significant principal components, [default: 20].")
    parser.add_argument('--clusterdim',metavar='INT',type=int,default=20,help='The principal components used for clustering, [default: 20].')
    parser.add_argument('--doubletpercentage',metavar='FLOAT',type=float,default=0.05,help='Assuming doublet formation rate, tailor for your dataset, [default: 0.05].')
    parser.add_argument('--mitpercentage',metavar='INT',type=int,default=15,help='Filter cells with mtgenes percentage, [default: 15].')
    parser.add_argument('--mtgenes',metavar='LIST',default='auto',help='Set mitochondrial genes(mtgene list file path) or false, [default: auto].')
    parser.add_argument('--minfeatures',metavar='INT',type=int,default=200,help='Filter cells with minimum nfeatures, [default: 200].',)
    parser.add_argument('--PCusage',metavar='INT',type=int,default=50,help='The total number of principal components for PCA, [default: 50].')
    parser.add_argument('--resolution',metavar='FLOAT',type=float,default=0.5,help='Cluster resolution, [default: 0.5].',)
    parser.add_argument('--species',metavar='STR',type=str,default='other',help='select species for cell annotation, only Human and Mouse can do auto annotation.')
    parser.add_argument('--outdir',metavar='DIR',help='output dir, [default: current directory].',default=os.getcwd())
    return parser