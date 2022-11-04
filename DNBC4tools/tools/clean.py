import os
from subprocess import check_call
from .utils import change_path,python_path
from DNBC4tools.__init__ import _root_dir

class Clean:
    def __init__(self, args):
        self.name = args.name
        self.outdir = args.outdir
        self.combine = args.combine

    def run(self):
        change_path()
        new_python = python_path()

        clean_cmd = ['%s %s/rna/clean_sample.py --outdir %s'%(new_python,_root_dir,self.outdir)]
        if self.name:
            clean_cmd += ['--name %s'%self.name]
        if self.combine:
            clean_cmd += ['--combine']
        clean_cmd = ' '.join(clean_cmd)
        #print(clean_cmd)
        check_call(clean_cmd,shell=True)

def clean(args):
    Clean(args).run()

def parse_clean(parser):
    parser.add_argument('--name',metavar='NAME',help='Sample name, [default is all sample in current directory]')
    parser.add_argument('--outdir',metavar='DIR',help='output dir, [default: current directory].', default=os.getcwd())
    parser.add_argument('--combine', action='store_true',help="Combine sample result.")

