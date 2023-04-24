import os
from dnbc4tools.tools.utils import rm_temp,str_mkdir,start_print_cmd,judgeFilexits,change_path,logging_call,read_json,bin_path
from dnbc4tools.__init__ import __root_dir__

class Data:
    def __init__(self, args):
        self.cDNAr1 = args.cDNAfastq1
        self.cDNAr2 = args.cDNAfastq2
        self.oligor1 = args.oligofastq1
        self.oligor2 = args.oligofastq2
        self.threads = args.threads
        self.name = args.name
        self.chemistry = args.chemistry
        self.darkreaction = args.darkreaction
        self.customize = args.customize
        self.outdir = os.path.join(args.outdir,args.name)
        self.genomeDir = args.genomeDir
        self.no_introns = args.no_introns

    def run(self):
        judgeFilexits(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.genomeDir)
        str_mkdir('%s/01.data'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        bin_command = bin_path()

        genomeDir = os.path.abspath(self.genomeDir)
        indexConfig = read_json('%s/ref.json'%genomeDir)
        gtf = indexConfig['gtf']
        judgeFilexits(gtf)
        chrmt = indexConfig['chrmt']

        # os.chdir('%s/log'%self.outdir) 
        scRNA_parse_cmd = ['%s/python %s/rna/src/star_anno.py --name %s --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --threads %s --chemistry %s --darkreaction %s --outdir %s --genomeDir %s --gtf %s '\
            %(bin_command,__root_dir__,self.name,self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.threads,self.chemistry,self.darkreaction,self.outdir,genomeDir,gtf)]
        if chrmt != 'None':
            scRNA_parse_cmd += ['--chrmt %s'%chrmt]
        if self.no_introns:
            scRNA_parse_cmd += ['--no_introns']
        if self.customize:
            scRNA_parse_cmd += ['--customize %s'%self.customize]
        scRNA_parse_cmd = ' '.join(scRNA_parse_cmd)
        start_print_cmd(scRNA_parse_cmd)
        final_sort_cmd = '%s/samtools sort -@ %s %s/01.data/final.bam -o %s/01.data/final_sorted.bam'\
            %(bin_command,self.threads,self.outdir,self.outdir)
        logging_call(final_sort_cmd,'data',self.outdir)
        str_mkdir('%s/01.data/raw_matrix'%self.outdir)
        raw_matrix_cmd = '%s/software/PISA count -@ %s -cb CB -anno-tag GN -umi UB -outdir %s/01.data/raw_matrix %s/01.data/final_sorted.bam'\
            %(__root_dir__,self.threads,self.outdir,self.outdir)
        print('\nGenerating raw expression matrix.')
        logging_call(raw_matrix_cmd,'data',self.outdir)
        rm_temp('%s/01.data/Aligned.out.bam'%self.outdir)
        rm_temp('%s/01.data/final.bam'%self.outdir)

def data(args):
    Data(args).run()

def helpInfo_data(parser):
    parser.add_argument(
        '--name', 
        metavar='NAME',
        help='Sample name.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output diretory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '--cDNAfastq1', 
        metavar='FASTQ',
        help='cDNA R1 fastq file, use commas to separate multiple files.', 
        required=True
        )
    parser.add_argument(
        '--cDNAfastq2', 
        metavar='FASTQ',
        help='cDNA R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with cDNAfastq1.', 
        required=True
        )
    parser.add_argument(
        '--oligofastq1', 
        metavar='FASTQ',
        help='oligo R1 fastq file, use commas to separate multiple files.',
        required=True
        )
    parser.add_argument(
        '--oligofastq2', 
        metavar='FASTQ',
        help='oligo R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with oligofastq1.',
        required=True
        )
    parser.add_argument(
        '--chemistry',
        metavar='STR',
        choices=["scRNAv1HT","scRNAv2HT","auto"],
        help='Chemistry version. Automatic detection is recommended. If setting, needs to be used with --darkreaction, can be "scRNAv1HT", "scRNAv2HT", [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--darkreaction',
        metavar='STR',
        help='Sequencing dark reaction. Automatic detection is recommended. If setting, needs to be used with --chemistry, use comma to separate cDNA and oligo, can be "R1,R1R2", "R1,R1", "unset,unset", [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--customize',
        metavar='STR',
        help='Customize files for whitelist and readstructure in JSON format for cDNA and oligo, use comma to separate cDNA and oligo.'
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads to use, [default: 4].'
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Not include intronic reads in count.'
        )
    return parser
