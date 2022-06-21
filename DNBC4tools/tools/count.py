import os
from .utils import str_mkdir,python_path,logging_call,judgeFilexits,change_path,rscript_path
from DNBC4tools.__init__ import _root_dir

class Count:
    def __init__(self,args):
        self.name = args.name
        self.bam = args.bam
        self.raw_matrix = args.raw_matrix
        self.cDNAbarcodeCount = args.cDNAbarcodeCount
        self.Indexreads = args.Indexreads
        self.oligobarcodeCount = args.oligobarcodeCount
        self.thread = args.thread
        self.oligotype = args.oligotype
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.minumi = args.minumi
        self.outdir = os.path.join(args.outdir,args.name)
    
    def run(self):
        judgeFilexits(self.bam,self.cDNAbarcodeCount,self.Indexreads,self.oligobarcodeCount,self.oligotype,self.raw_matrix)
        str_mkdir('%s/02.count'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        new_python = python_path()
        new_rscript = rscript_path()

        cellCalling_cmd = '%s %s/rna/cell_calling.R --matrix %s --outdir %s/02.count/ --method %s --expectcells %s --forcecells %s --minumi %s'\
            %(new_rscript,_root_dir,self.raw_matrix,self.outdir,self.calling_method,self.expectcells,self.forcecells,self.minumi)
        get_barcode_cmd = '%s %s/rna/get_barcode.py --raw %s --select %s/02.count/beads_barcodes_hex.txt -o %s/02.count/'\
            %(new_python,_root_dir,self.cDNAbarcodeCount,self.outdir,self.outdir)
        mergeBarcodes_cmd = '%s/soft/mergeBarcodes -b %s/02.count/beads_barcode_all.txt -f %s -n %s -o %s/02.count/'\
            %(_root_dir,self.outdir,self.Indexreads,self.name,self.outdir)
        similiarBeads_cmd = '%s/soft/s1.get.similarityOfBeads -n %s %s %s/02.count/%s_CB_UB_count.txt %s/02.count/beads_barcodes.txt %s %s/02.count/Similarity.all.csv %s/02.count/Similarity.droplet.csv %s/02.count/Similarity.droplet.filtered.csv'\
            %(_root_dir,self.thread,self.name,self.outdir,self.name,self.outdir,self.oligotype,self.outdir,self.outdir,self.outdir)
        combineBeads_cmd = '%s %s/rna/combinedListOfBeads.py --similarity_droplet %s/02.count/Similarity.droplet.csv --beads_list %s/02.count/beads_barcodes.txt --combined_list %s/02.count/%s_combined_list.txt'\
            %(new_python,_root_dir,self.outdir,self.outdir,self.outdir,self.name)
        CellMerge_cmd = '%s %s/rna/cellMerge.py --indir %s/02.count --name %s'\
            %(new_python,_root_dir,self.outdir,self.name)
        tagAdd_cmd = '%s/soft/tagAdd -n %s -bam %s -file %s/02.count/%s_barcodeTranslate_hex.txt -out %s/02.count/anno_decon.bam -tag_check CB:Z: -tag_add DB:Z: -root %s'\
            %(_root_dir,self.thread,self.bam,self.outdir,self.name,self.outdir,_root_dir)
        PISA_count_cmd = '%s/soft/PISA count -one-hit -@ %s -cb DB -anno-tag GN -umi UB -outdir %s/02.count/filter_matrix %s/02.count/anno_decon.bam'\
            %(_root_dir,self.thread,self.outdir,self.outdir)
        cellReport_cmd = '%s %s/rna/cell_report.R -M %s/02.count/filter_matrix -O %s/02.count'\
            %(new_rscript,_root_dir,self.outdir,self.outdir)
        bam_sort_cmd = 'samtools sort -@ %s %s/02.count/anno_decon.bam -o %s/02.count/anno_decon_sorted.bam'\
            %(self.thread,self.outdir,self.outdir)
        bam_index_cmd = 'samtools index -@ %s %s/02.count/anno_decon_sorted.bam'\
            %(self.thread,self.outdir)
        saturation_cmd = '%s %s/rna/saturation.py -i %s/02.count/anno_decon_sorted.bam -o %s/02.count -f %s/02.count/cellCount_report.csv --quality 20 --threads %s'\
            %(new_python,_root_dir,self.outdir,self.outdir,self.outdir,self.thread)
        

        logging_call(cellCalling_cmd,'count',self.outdir)
        logging_call(get_barcode_cmd,'count',self.outdir)
        logging_call(mergeBarcodes_cmd,'count',self.outdir)
        logging_call(similiarBeads_cmd,'count',self.outdir)
        logging_call(combineBeads_cmd,'count',self.outdir)
        logging_call(CellMerge_cmd,'count',self.outdir)
        logging_call(tagAdd_cmd,'count',self.outdir)
        str_mkdir('%s/02.count/filter_matrix'%self.outdir)
        logging_call(PISA_count_cmd,'count',self.outdir)
        logging_call(cellReport_cmd,'count',self.outdir)
        logging_call(bam_sort_cmd,'count',self.outdir)
        logging_call(bam_index_cmd,'count',self.outdir)
        logging_call(saturation_cmd,'count',self.outdir)
        
        #rm_temp('%s/02.count/%s_CB_UB_count.txt'%(self.outdir,self.name),'%s/02.count/beads_barcode_all.txt'%self.outdir)

def count(args):
    Count(args).run()

def parse_count(parser):
    parser.add_argument('--name',metavar='NAME',help='sample name.')
    parser.add_argument('--bam',metavar='FILE',help='Bam file after star and anno, eg./01.data/final.bam.')
    parser.add_argument('--raw_matrix',metavar='FILE',help='Raw matrix dir, eg./01.data/raw_matrix.')
    parser.add_argument('--cDNAbarcodeCount',metavar='FILE',help='Read count per cell barcode for cDNA, eg./01.data/cDNA_barcode_counts_raw.txt.',)
    parser.add_argument('--Indexreads',metavar='FILE',help='Barcode reads generate by scRNAparse, eg./01.data/Index_reads.fq.gz.')
    parser.add_argument('--oligobarcodeCount',metavar='FILE',help='Read count per cell barcode for oligo, eg./01.data/Index_barcode_counts_raw.txt.')
    parser.add_argument('--oligotype',metavar='FILE',help='Whitelist for oligo index, [default: oligo_type8.txt].',default='%s/config/oligo_type8.txt'%_root_dir)
    parser.add_argument('--thread',metavar='INT',help='Analysis threads. [default: 4].',type=int,default=4)
    parser.add_argument('--outdir',metavar='DIR',help='output dir, [default: current directory].',default=os.getcwd())
    parser.add_argument('--calling_method',metavar='STR',help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].', default='emptydrops')
    parser.add_argument('--expectcells',metavar='INT',help='Expected number of recovered beads, used as input to cell calling algorithm, [default: 3000].', default=3000)
    parser.add_argument('--forcecells',metavar='INT',help='Force pipeline to use this number of beads, bypassing cell calling algorithm.',default=0)
    parser.add_argument('--minumi',metavar='INT',help='The min umi for use emptydrops, [default: 1000].', default=1000)
    return parser