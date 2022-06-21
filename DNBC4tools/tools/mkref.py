import os
from .utils import str_mkdir,start_print_cmd,judgeFilexits,change_path,python_path
from DNBC4tools.__init__ import _root_dir

class Ref:
    def __init__(self, args):
        self.action = args.action
        self.ingtf = args.ingtf
        self.outgtf = args.outgtf
        self.attribute = args.attribute 
        self.outstat = args.outstat
        self.type = args.type
        self.fasta = args.fasta
        self.star_dir = args.star_dir
        self.thread = args.thread

    def run(self):
        
        change_path()
        new_python = python_path()
        if self.action == 'mkgtf':
            judgeFilexits(self.ingtf)
            list = (" ".join(str(i) for i in self.attribute))
            mkgtf_cmd = '%s %s/rna/star_ref.py --action mkgtf --ingtf %s --outgtf %s --attribute %s'\
                %(new_python,_root_dir,self.ingtf,self.outgtf,list)
            start_print_cmd(mkgtf_cmd)
        if self.action == 'stat':
            judgeFilexits(self.ingtf)
            stat_cmd = '%s %s/rna/star_ref.py --action stat --ingtf %s --type %s --outstat %s'\
                %(new_python,_root_dir,self.ingtf,self.type,self.outstat)
            start_print_cmd(stat_cmd)
        if self.action == 'mkref':
            judgeFilexits(self.ingtf,self.fasta)
            str_mkdir('%s'%self.star_dir)
            mkref_cmd = '%s %s/rna/star_ref.py --action mkref --ingtf %s --fasta %s --star_dir %s --star %s/soft/scStar --thread %s'\
                %(new_python,_root_dir,self.ingtf,self.fasta,self.star_dir,_root_dir,self.thread)
            start_print_cmd(mkref_cmd)

def mkref(args):
    Ref(args).run()

def parse_mkref(parser):
    parser.add_argument('--action',metavar='SELECT',default='mkref', choices=['mkref', 'mkgtf', 'stat'], help='Select the action for your program, include mkref,mkgtf,stat, [default: mkref].')
    parser.add_argument('--ingtf', metavar='FILE' ,help='Set ingtf in mkref,mkgtf or stat.')
    parser.add_argument('--outgtf',metavar='FILE', help='Set outgtf in mkgtf.')
    parser.add_argument('--attribute',metavar='DICT',default = ['gene_type:protein_coding'], nargs='+',help='Set the filter parameter in mkgtf, Key-value pair in attributes field to be kept in the GTF, default is gene_type:protein_coding, you can set up multiple groups separated by blankspace.')
    parser.add_argument('--outstat',metavar='FILE', default = 'gtf_type.txt',help='Set the stats outfile in stat, [default: "gtf_type.txt" in current dir].')
    parser.add_argument('--type',metavar='STR', default = 'gene_type',help='Set the the type for stat, [default: gene_type].')
    parser.add_argument('--fasta',metavar='FASTA',help='Set the fasta in mkref.')
    parser.add_argument('--star_dir',metavar='DIR',default=os.getcwd(),help='Set the star indexdir in mkref, [default: current dir].')
    parser.add_argument('--thread',metavar='INT', default=4,help='Set the threads in mkref, [default: 4].')
    return parser