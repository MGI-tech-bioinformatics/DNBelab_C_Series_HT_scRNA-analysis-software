import os,shutil
from dnbc4tools.tools.utils import str_mkdir,logging_call,change_path,read_json,judgeFilexits,rm_temp
from dnbc4tools.__init__ import __root_dir__

class Report:
    def __init__(self,args):
        self.name = args.name
        self.genomeDir = args.genomeDir
        self.threads = args.threads
        self.no_introns = args.no_introns
        self.outdir = os.path.join(args.outdir,args.name)

    def run(self):
        rm_temp('%s/04.report/div/anno.div'%self.outdir,'%s/04.report/div/anno_chsize.div'%self.outdir,'%s/04.report/div/anno.html'%self.outdir)
        str_mkdir('%s/04.report/div'%self.outdir)
        str_mkdir('%s/04.report/base64'%self.outdir)
        str_mkdir('%s/04.report/table'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        genomeDir = os.path.abspath(self.genomeDir)
        indexConfig = read_json('%s/ref.json'%genomeDir)
        species = indexConfig['species']

        ### plotly
        print('\nSummarizing results')
        from dnbc4tools.tools._plotly import plot_jaccard_knee_frag,plot_cluster,plot_cluster_umi,plot_clusteranno,plot_saturation
        judgeFilexits(
            '%s/02.count/cutoff.csv'%self.outdir,
            '%s/02.count/saturation.xls'%self.outdir,
            '%s/03.analysis/cluster.csv'%self.outdir
        )
        plot_jaccard_knee_frag('%s/02.count/cutoff.csv'%self.outdir,'%s/04.report/div'%self.outdir)
        plot_cluster('%s/03.analysis/cluster.csv'%self.outdir,'%s/04.report/div'%self.outdir)
        plot_cluster_umi('%s/03.analysis/cluster.csv'%self.outdir,'%s/04.report/div'%self.outdir)
        with open('%s/03.analysis/cluster.csv'%self.outdir,'r') as clusterfile:
            firstline = clusterfile.readline().rstrip()
        if 'Predicted cell type' in firstline:
            plot_clusteranno('%s/03.analysis/cluster.csv'%self.outdir,'%s/04.report/div'%self.outdir)
        plot_saturation('%s/02.count/saturation.xls'%self.outdir,'%s/04.report/div'%self.outdir)

        ### png2base
        from dnbc4tools.tools.utils import png_to_base64,csv_datatable
        judgeFilexits(
            '%s/02.count/cellNumber_merge.png'%self.outdir,
            '%s/03.analysis/raw_QCplot.png'%self.outdir,
            '%s/03.analysis/marker.csv'%self.outdir,
        )
        pictures = {'%s/02.count/cellNumber_merge.png'%self.outdir:'6','%s/03.analysis/raw_QCplot.png'%self.outdir:'7'}
        for infile,outfile in pictures.items():
            png_to_base64(infile,outfile,'%s/04.report/base64'%self.outdir)

        ### marker table
        csv_datatable('%s/03.analysis/marker.csv'%self.outdir,'%s/04.report/table/marker-table.txt'%self.outdir)

        ### html 
        intron = 'False' if self.no_introns else 'True' 
        from dnbc4tools.rna.src.generate_report import write_param_to_template
        report_file, metrics_summary_df = write_param_to_template(
            '%s/template/template_scRNA.html'%__root_dir__,
            self.name,
            self.outdir,
            intron,
            species
            )
        htmlFile = open('%s/04.report/%s_scRNA_report.html'%(self.outdir,self.name),'w')
        htmlFile.write(report_file)
        htmlFile.close()
        metrics_summary_df.to_csv('%s/04.report/metrics_summary.xls'%self.outdir,sep='\t',index=None)

        ### copy file
        str_mkdir('%s/output'%self.outdir)
        if os.path.exists('%s/04.report/metrics_summary.xls'%self.outdir):
            shutil.copy("%s/04.report/metrics_summary.xls"%self.outdir,'%s/output'%self.outdir)
        if os.path.exists('%s/04.report/%s_scRNA_report.html'%(self.outdir,self.name)):
            shutil.copy("%s/04.report/%s_scRNA_report.html"%(self.outdir,self.name),'%s/output'%self.outdir)
        if os.path.exists('%s/02.count/filter_feature.h5ad'%self.outdir):
            shutil.copy('%s/02.count/filter_feature.h5ad'%self.outdir,'%s/output'%self.outdir)
        if os.path.exists("%s/02.count/filter_matrix"%self.outdir):
            if os.path.exists('%s/output/filter_matrix'%self.outdir):
                shutil.rmtree('%s/output/filter_matrix'%self.outdir)
            shutil.copytree("%s/02.count/filter_matrix"%self.outdir,'%s/output/filter_matrix'%self.outdir,dirs_exist_ok=True)
        if os.path.exists("%s/01.data/raw_matrix"%self.outdir):
            if os.path.exists('%s/output/raw_matrix'%self.outdir):
                shutil.rmtree('%s/output/raw_matrix'%self.outdir)
            shutil.copytree("%s/01.data/raw_matrix"%self.outdir,'%s/output/raw_matrix'%self.outdir,dirs_exist_ok=True)
        
        ### move bam
        if os.path.exists('%s/02.count/anno_decon_sorted.bam'%self.outdir):
            if(os.path.exists("%s/output/anno_decon_sorted.bam"%self.outdir)):
                os.remove("%s/output/anno_decon_sorted.bam"%self.outdir)
            shutil.move("%s/02.count/anno_decon_sorted.bam"%self.outdir,'%s/output'%self.outdir)
        if os.path.exists('%s/02.count/anno_decon_sorted.bam.bai'%self.outdir):
            if(os.path.exists("%s/output/anno_decon_sorted.bam.bai"%self.outdir)):
                os.remove("%s/output/anno_decon_sorted.bam.bai"%self.outdir)
            shutil.move("%s/02.count/anno_decon_sorted.bam.bai"%self.outdir,'%s/output'%self.outdir)
        if os.path.exists('%s/02.count/anno_decon_sorted.bam.csi'%self.outdir):
            if(os.path.exists("%s/output/anno_decon_sorted.bam.csi"%self.outdir)):
                os.remove("%s/output/anno_decon_sorted.bam.csi"%self.outdir)
            shutil.move("%s/02.count/anno_decon_sorted.bam.csi"%self.outdir,'%s/output'%self.outdir)

        ### other matrix
        splice_matrix_cmd = '%s/software/PISA count -one-hit -@ %s -cb DB -ttype E,S -anno-tag GN -umi UB -list %s/02.count/cell.id -outdir %s/output/attachment/splice_matrix %s/output/anno_decon_sorted.bam' \
            %(__root_dir__,self.threads,self.outdir,self.outdir,self.outdir)
        RNAvelocity_matrix_cmd = '%s/software/PISA count -one-hit -@ %s -cb DB -velo -anno-tag GN -umi UB -list %s/02.count/cell.id -outdir %s/output/attachment/RNAvelocity_matrix %s/output/anno_decon_sorted.bam' \
            %(__root_dir__,self.threads,self.outdir,self.outdir,self.outdir)
        if not self.no_introns:
            str_mkdir('%s/output/attachment/splice_matrix'%self.outdir)
            str_mkdir('%s/output/attachment/RNAvelocity_matrix'%self.outdir)
            logging_call(splice_matrix_cmd,'report',self.outdir)
            logging_call(RNAvelocity_matrix_cmd,'report',self.outdir)
        else:
            if os.path.exists('%s/output/attachment/splice_matrix'%self.outdir):
                shutil.rmtree('%s/output/attachment/splice_matrix'%self.outdir)
            if os.path.exists('%s/output/attachment/RNAvelocity_matrix'%self.outdir):
                shutil.rmtree('%s/output/attachment/RNAvelocity_matrix'%self.outdir)

def report(args):
    Report(args).run()

def helpInfo_report(parser):
    parser.add_argument(
        '--name',
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
        help='output dir, [default: current directory].',
        default=os.getcwd()
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Analysis threads, [default: 4].'
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Not include intronic reads in count.'
        )
    return parser
