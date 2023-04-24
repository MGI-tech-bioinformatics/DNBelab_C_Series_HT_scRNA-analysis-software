import os
from dnbc4tools.tools.utils import str_mkdir,logging_call,judgeFilexits,change_path,bin_path
from dnbc4tools.__init__ import __root_dir__

def get_barcode(raw,hex_barcode):
    from dnbc4tools.tools.utils import seq_comp
    import pandas as pd
    barcode_all = pd.read_table(raw,sep = '\t',header=None)
    barcode_all.columns = ['barcode','count']
    barcode_all["hex"]= barcode_all["barcode"].map(seq_comp)
    
    select_barcode = []
    with open(hex_barcode,'r') as select:
        for line in select:
            line = line.strip()
            select_barcode.append(line)
    select_df = barcode_all.loc[barcode_all['hex'].isin(select_barcode)]
    
    return barcode_all,select_df

def matrix_summary(matrixpath,outdir,cellreport):
    from dnbc4tools.tools.utils import read_anndata
    import scanpy as sc
    adata = read_anndata(matrixpath)
    adata.write("%s/filter_feature.h5ad"%outdir)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    total_gene = str(adata.var.shape[0])
    mean_gene = str(round(adata.obs['n_genes_by_counts'].mean()))
    median_gene = str(round(adata.obs['n_genes_by_counts'].median()))
    mean_umi = str(round(adata.obs['total_counts'].mean()))
    median_umi = str(round(adata.obs['total_counts'].median()))
    with open(cellreport,'a') as reportfile:
        reportfile.write('Mean UMI counts per cell,%s'%mean_umi+'\n')
        reportfile.write('Median UMI Counts per Cell,%s'%median_umi+'\n')
        reportfile.write('Total Genes Detected,%s'%total_gene+'\n')
        reportfile.write('Mean Genes per Cell,%s'%mean_gene+'\n')
        reportfile.write('Median Genes per Cell,%s'%median_gene+'\n')

class Count:
    def __init__(self,args):
        self.name = args.name
        self.threads = args.threads
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.minumi = args.minumi
        self.outdir = os.path.join(args.outdir,args.name)
    
    def run(self):
        judgeFilexits(
            '%s/01.data/final_sorted.bam'%self.outdir,
            '%s/01.data/cDNA_barcode_counts_raw.txt'%self.outdir,
            '%s/01.data/Index_reads.fq.gz'%self.outdir,
            '%s/01.data/beads_stat.txt'%self.outdir,
            '%s/01.data/raw_matrix'%self.outdir
            )
        str_mkdir('%s/02.count'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        bin_command = bin_path()

        ## cell calling using DropletUtils
        print('\nCalling cell barcodes.')
        cellCalling_cmd = '%s/Rscript %s/rna/src/cell_calling.R --matrix %s/01.data/raw_matrix --outdir %s/02.count/ --method %s --expectcells %s --forcecells %s --minumi %s'\
            %(bin_command,__root_dir__,self.outdir,self.outdir,self.calling_method,self.expectcells,self.forcecells,self.minumi)
        logging_call(cellCalling_cmd,'count',self.outdir)

        ### get all barcode and select barcode
        barcode_all,select_df = get_barcode('%s/01.data/cDNA_barcode_counts_raw.txt'%self.outdir,
                                            '%s/02.count/beads_barcodes_hex.txt'%self.outdir)
        barcode_all['barcode'].to_csv(os.path.join(self.outdir,'02.count/beads_barcode_all.txt'),index=False,header=False)
        select_df['barcode'].to_csv(os.path.join(self.outdir,'02.count/beads_barcodes.txt'),index=False,header=False)

        print('\nCalculating bead similarity and merging beads..')
        ### using index reads to merge beads
        mergeBarcodes_cmd = '%s/software/mergeBarcodes -b %s/02.count/beads_barcode_all.txt -f %s/01.data/Index_reads.fq.gz -n %s -o %s/02.count/'\
            %(__root_dir__,self.outdir,self.outdir,self.name,self.outdir)
        similiarBeads_cmd = '%s/software/similarityOfBeads -n %s %s %s/02.count/%s_CB_UB_count.txt %s/02.count/beads_barcodes.txt %s/config/oligo_type8.txt %s/02.count/Similarity.all.csv %s/02.count/Similarity.droplet.csv %s/02.count/Similarity.droplet.filtered.csv'\
            %(__root_dir__,self.threads,self.name,self.outdir,self.name,self.outdir,__root_dir__,self.outdir,self.outdir,self.outdir)
        logging_call(mergeBarcodes_cmd,'count',self.outdir)
        logging_call(similiarBeads_cmd,'count',self.outdir)

        ### merge beads list
        from dnbc4tools.rna.src.combinedListOfBeads import similarity_droplet_file
        similarity_droplet_file('%s/02.count/Similarity.droplet.csv'%self.outdir,
                                '%s/02.count/beads_barcodes.txt'%self.outdir,
                                '%s/02.count/%s_combined_list.txt'%(self.outdir,self.name))


        ### summary beads merge 
        from dnbc4tools.rna.src.cellMerge import summary_count
        summary_count('%s/02.count/%s_combined_list.txt'%(self.outdir,self.name),
                      '%s/02.count/beads_barcodes.txt'%self.outdir,
                      '%s/01.data/beads_stat.txt'%self.outdir,
                      '%s/02.count/%s_barcodeTranslate.txt'%(self.outdir,self.name),
                      '%s/02.count/%s_barcodeTranslate_hex.txt'%(self.outdir,self.name),
                      '%s/02.count/cell.id'%self.outdir,
                      '%s/02.count/cellCount_report.csv'%self.outdir,
                      '%s/02.count'%self.outdir)

        ### add DB tag for bam
        tagAdd_cmd = '%s/software/tagAdd -n %s -bam %s/01.data/final_sorted.bam -file %s/02.count/%s_barcodeTranslate_hex.txt -out %s/02.count/anno_decon_sorted.bam -tag_check CB:Z: -tag_add DB:Z: '\
            %(__root_dir__,self.threads,self.outdir,self.outdir,self.name,self.outdir)
        logging_call(tagAdd_cmd,'count',self.outdir)

        ### PISA count matrix
        print('\nGenerating filter expression matrix.')
        str_mkdir('%s/02.count/filter_matrix'%self.outdir)
        PISA_count_cmd = '%s/software/PISA count -one-hit -@ %s -cb DB -anno-tag GN -umi UB -list %s/02.count/cell.id -outdir %s/02.count/filter_matrix %s/02.count/anno_decon_sorted.bam'\
            %(__root_dir__,self.threads,self.outdir,self.outdir,self.outdir)
        logging_call(PISA_count_cmd,'count',self.outdir)

        ### get cell report
        matrix_summary('%s/02.count/filter_matrix'%self.outdir,
                        '%s/02.count'%self.outdir,
                        '%s/02.count/cellCount_report.csv'%self.outdir)

        ### get bam index
        def create_index(threads,bam):
            try:
                bam_index_cmd = '%s/samtools index -@ %s %s'%(bin_command,threads,bam)
                logging_call(bam_index_cmd,'count',self.outdir)
            except Exception as e:
                print('build csi index for bam')
                bam_index_cmd = '%s/samtools index -c -@ %s %s'%(bin_command,threads,bam)
                logging_call(bam_index_cmd,'count',self.outdir)
        create_index(self.threads,'%s/02.count/anno_decon_sorted.bam'%self.outdir)

        ### cal saturation
        saturation_cmd = '%s/python %s/rna/src/saturation.py -i %s/02.count/anno_decon_sorted.bam -o %s/02.count -f %s/02.count/cellCount_report.csv --quality 20 --threads %s'\
            %(bin_command,__root_dir__,self.outdir,self.outdir,self.outdir,self.threads)
        logging_call(saturation_cmd,'count',self.outdir)


def count(args):
    Count(args).run()

def helpInfo_count(parser):
    parser.add_argument(
        '--name',
        metavar='NAME',
        help='sample name.'
        )
    parser.add_argument(
        '--threads',
        metavar='INT',
        help='Analysis threads. [default: 4].',
        type=int,default=4
        )
    parser.add_argument(
        '--outdir',
        metavar='DIR',
        help='output dir, [default: current directory].',
        default=os.getcwd()
        )
    parser.add_argument(
        '--calling_method',
        metavar='STR',
        help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].', 
        default='emptydrops'
        )
    parser.add_argument(
        '--expectcells',
        metavar='INT',
        help='Expected number of recovered beads, used as input to cell calling algorithm, [default: 3000].', 
        default=3000
        )
    parser.add_argument(
        '--forcecells',
        metavar='INT',
        help='Force pipeline to use this number of beads, bypassing cell calling algorithm.',
        default=0
        )
    parser.add_argument(
        '--minumi',
        metavar='INT',
        help='The min umi for use emptydrops, [default: 1000].', 
        default=1000
        )
    return parser