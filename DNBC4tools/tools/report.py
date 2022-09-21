import os,shutil
from .utils import str_mkdir,logging_call,change_path,python_path
from DNBC4tools.__init__ import _root_dir

class Report:
    def __init__(self,args):
        self.name = args.name
        self.species = args.species
        self.thread = args.thread
        self.no_bam = args.no_bam
        self.no_introns = args.no_introns
        self.outdir = os.path.join(args.outdir,args.name)
    def run(self):
        str_mkdir('%s/04.report'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        new_python = python_path()
        pre_cmd = '%s %s/rna/pre_process.py --outPath %s --sample %s'\
            %(new_python,_root_dir,self.outdir,self.name)
        generate_report_cmd = ['%s %s/rna/generate_report.py --outPath %s --htmlTemplate %s/template/template.html --name %s --species %s' \
            %(new_python,_root_dir,self.outdir,_root_dir,self.name,self.species)]
        if self.no_introns:
            generate_report_cmd += ['--intron false']
        else:
            generate_report_cmd += ['--intron true']
        generate_report_cmd = ' '.join(generate_report_cmd)
        generate_output_cmd = ['%s %s/rna/report_output.py --indir %s'%(new_python,_root_dir,self.outdir)]
        if self.no_bam:
            generate_output_cmd += ['--no_bam']
        generate_output_cmd = ' '.join(generate_output_cmd)

        if self.no_bam:
            anno_bam_dir = '02.count'
        else:
            anno_bam_dir = 'output'
        splice_matrix_cmd = '%s/soft/PISA count -one-hit -@ %s -cb DB -ttype E,S -anno-tag GN -umi UB -list %s/02.count/cell.id -outdir %s/output/attachment/splice_matrix %s/%s/anno_decon_sorted.bam' \
            %(_root_dir,self.thread,self.outdir,self.outdir,self.outdir,anno_bam_dir)
        RNAvelocity_matrix_cmd = '%s/soft/PISA count -one-hit -@ %s -cb DB -velo -anno-tag GN -umi UB -list %s/02.count/cell.id -outdir %s/output/attachment/RNAvelocity_matrix %s/%s/anno_decon_sorted.bam' \
            %(_root_dir,self.thread,self.outdir,self.outdir,self.outdir,anno_bam_dir)
        
        logging_call(pre_cmd,'report',self.outdir)
        logging_call(generate_report_cmd,'report',self.outdir)
        logging_call(generate_output_cmd,'report',self.outdir)
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

def parse_report(parser):
    parser.add_argument('--name',required=True,help='Sample name.')
    parser.add_argument('--species',type=str,default='NA',help='Species name.')
    parser.add_argument('--outdir',help='output dir, [default: current directory].',default=os.getcwd())
    parser.add_argument('--thread',type=int, metavar='INT',default=4,help='Analysis threads, [default: 4].')
    parser.add_argument('--no_bam', action='store_true',help='Do not move filter bam file to output dir.')
    parser.add_argument('--no_introns', action='store_true',help='Not include intronic reads in count.')
    return parser
