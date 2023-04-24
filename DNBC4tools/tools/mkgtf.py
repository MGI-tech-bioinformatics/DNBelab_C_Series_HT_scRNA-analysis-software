import re
import collections
from typing import List, Iterable

from collections import Counter
from dnbc4tools.tools.utils import judgeFilexits, change_path
from dnbc4tools.__init__ import __root_dir__


def read_gtf(fp: str, feature: str) -> Iterable[List[str]]:
    '''
    Reads the gtf file and returns an iterable of lists.
    Each list contains lines of a single gene.
    '''
    lines = []
    with open(fp) as f:
        for line in f:
            newline = line.strip()
            if newline.startswith('#'):
                lines.append(line)
            elif newline == '':
                continue
            else:
                lst = newline.split('\t')
                if lst[2] == feature:
                    yield lines
                    lines = []
                    lines.append(line)
                else:
                    lines.append(line)
        yield lines


def filtergtf(gtf: str, filtergtf: str, keyword: str, attribute: List[str]):
    '''
    Filters gtf file using a keyword and a list of attributes.
    '''
    d = dict([(key,keyword) for key in attribute])
    gtfread = read_gtf(gtf,'gene')
    result = open(filtergtf,'w')
    for i in gtfread:
        if i:
            if i[0].startswith('#'):
                result.writelines(i)
            else:
                aDict = collections.OrderedDict()
                pattern = re.compile(r'(\S+?)\s*"(.*?)"')
                for m in re.finditer(pattern, i[0].split('\t')[-1]):
                    key = m.group(1)
                    value = m.group(2)
                    aDict[key] = value
                for key1,value1 in aDict.items():
                    for key2,value2 in d.items():
                        if key1 == value2 and key2 == value1:
                            result.writelines(i)
    result.close()

def statgtf(gtf: str,keyword: str,outfile: str):
    '''
    Count the occurrence of the given keyword in the attributes of the GTF file, 
    and output the results to the outfile
    '''
    with open(gtf,'r') as fp:
        sumDict = []
        for line in fp:
            line = line.strip()
            if line.startswith("#"):
                continue
            elif line == '':
                continue
            else:
                lst = line.split('\t')
                if lst[2] == 'gene':
                    aDict = collections.OrderedDict()
                    pattern = re.compile(r'(\S+?)\s*"(.*?)"')
                    for m in re.finditer(pattern, lst[-1]):
                        key = m.group(1)
                        value = m.group(2)
                        aDict[key] = value
                    sumDict.append(aDict[keyword])
        result = Counter(sumDict)
        outfile = open(outfile,'w')
        outfile.write('Type'+'\t'+'Count'+'\n')
        for k,v in sorted(result.items(), key = lambda x:x[1], reverse=True):
            outfile.write(f'{k}\t{v}\n')
        outfile.close()


normal_include = 'protein_coding,\
lncRNA,\
lincRNA,\
antisense,\
IG_C_gene,\
IG_D_gene,\
IG_J_gene,\
IG_LV_gene,\
IG_V_gene,\
IG_V_pseudogene,\
IG_J_pseudogene,\
IG_C_pseudogene,\
TR_C_gene,\
TR_D_gene,\
TR_J_gene,\
TR_V_gene,\
TR_V_pseudogene,\
TR_J_pseudogene'


class Mkgtf:
    def __init__(self, args):
        """
        Constructor for Mkgtf class.
        Args:
        args (argparse.Namespace): Namespace object that contains parsed command-line arguments.
        """
        self.action: str = args.action
        self.ingtf: str = args.ingtf
        self.output: str = args.output
        self.include: str = args.include
        self.type: str = args.type
        
    def run(self):
        """
        Run the specified action according to the parsed command-line arguments.
        """
        change_path()

        if self.action == 'mkgtf':
            judgeFilexits(self.ingtf)
            gene_types: List[str] = [i.strip() for i in self.include.split(',')]
            print(f"\033[0;32;40mUsing \"{self.type}\" for Analysis\033[0m")
            print("\033[0;32;40mThe types of genes obtained are:\033[0m")
            for gene_type in gene_types:
                print(gene_type)
            filtergtf(self.ingtf, self.output, self.type, gene_types)
            print("\033[0;32;40mAnalysis Complete\033[0m")

        if self.action == 'stat':
            judgeFilexits(self.ingtf)
            print(f"\033[0;32;40mUsing \"{self.type}\" for Analysis\033[0m")
            statgtf(self.ingtf, self.type, self.output)
            print("\033[0;32;40mAnalysis Complete\033[0m")

def mkgtf(args):
    Mkgtf(args).run()

def helpInfo_mkgtf(parser):
    parser.add_argument(
        '--action',
        metavar='SELECT',
        default='mkgtf', 
        choices=['mkgtf', 'stat'], 
        help="Select the action for your program, including 'mkgtf' and 'stat', [default: mkgtf]."
        )
    parser.add_argument(
        '--ingtf', 
        metavar='FILE' ,
        help='Path to the GTF file with annotations.'
        )
    parser.add_argument(
        '--output',
        metavar='FILE', 
        help='Path to output file.'
        )
    parser.add_argument(
        '--include',
        metavar='STR',
        default = normal_include,
        help="Set the filter parameter in 'mkgtf'. You can set up multiple filters separated by commas."
        )
    parser.add_argument(
        '--type',
        metavar='STR', 
        default = 'gene_biotype',
        help='Set according to the tag of gene type in GTF attributes, [default: gene_biotype].'
        )
    return parser
