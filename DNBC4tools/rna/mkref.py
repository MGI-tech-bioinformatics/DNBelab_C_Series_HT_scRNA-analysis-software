import os,math,re,json
from dnbc4tools.tools.utils import judgeFilexits,change_path,start_print_cmd
from dnbc4tools.__init__ import __root_dir__
from dnbc4tools.tools.mkgtf import *

def star_index(fasta,gtf,genomeDir,star_program,limitram,threads):
    if not os.path.exists(genomeDir):
        os.system('mkdir -p %s'%genomeDir)
    genome_size = os.path.getsize(fasta)
    SAindexNbases = min(int(math.log2(genome_size)/2 - 1),14)
    star_cmd = '%s --runMode genomeGenerate --runThreadN %s --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang 99 --limitGenomeGenerateRAM %s --genomeSAindexNbases %s '\
        %(star_program,threads,genomeDir,fasta,gtf,limitram,SAindexNbases)
    print('STAR verison: 2.7.2b')
    print('runMode: genomeGenerate')
    print('runThreadN: %s'%threads)
    print('limitGenomeGenerateRAM: %s'%limitram)
    print('genomeSAindexNbases: %s'%SAindexNbases)
    print('genomeDir: %s'%genomeDir)
    print('fasta: %s'%fasta)
    print('gtf: %s'%gtf)
    start_print_cmd(star_cmd)

def mtgenes_list(chrM,genomeDir,gtf):
    if chrM == 'auto':
        chrmtlist = ['chrM','MT','chrMT','mt','Mt']
        if os.path.exists(os.path.join(genomeDir,'chrName.txt')):
            with open(os.path.join(genomeDir,'chrName.txt'),'r') as file:
                chrlist =  file.read().splitlines()
                union = list(set(chrmtlist) & set(chrlist))
                if union:
                    chrMT = union[0]
                else:
                    chrMT = 'None'
        else:
            chrMT = 'None'
    else:
        chrMT = chrM

    if chrMT != 'None':
        MTgene = []
        with open(gtf,'r') as gtffile:
            for line in gtffile:
                line = line.strip()
                if line.startswith('#'):
                    pass
                else:
                    lst = line.split('\t')
                    if lst[0] == chrMT:
                        if lst[2] == "gene":
                            aDict = collections.OrderedDict()
                            pattern = re.compile(r'(\S+?)\s*"(.*?)"')
                            for m in re.finditer(pattern, lst[-1]):
                                key = m.group(1)
                                value = m.group(2)
                                aDict[key] = value
                            if "gene_name" in aDict:
                                MTgene.append(aDict['gene_name'])
                            else:
                                MTgene.append(aDict['gene_id'])
        if MTgene:
            result = open(os.path.join(genomeDir,'mtgene.list'),'w')
            for line in MTgene:
                result.write(line+'\n')
            result.close()
    return chrMT

def write_config(genomeDir,species,chrM,fasta,gtf):
    ref_dict = collections.OrderedDict()
    ref_dict['species'] = str(species)
    ref_dict['genome'] = os.path.abspath(fasta)
    ref_dict['gtf'] = os.path.abspath(gtf)
    ref_dict['genomeDir'] = os.path.abspath(genomeDir)
    chrMT = mtgenes_list(chrM,genomeDir,gtf)
    ref_dict['chrmt'] = str(chrMT)
    if os.path.exists(os.path.join(genomeDir,'mtgene.list')):
        ref_dict['mtgenes'] = os.path.abspath("%s/mtgene.list"%genomeDir)
    else:
        ref_dict['mtgenes'] = 'None'
    with open("%s/ref.json"%genomeDir, "w", encoding='utf-8') as jsonfile:
        json.dump(ref_dict, jsonfile, indent=4,  ensure_ascii=False) 
        jsonfile.write("\n")

class Ref:
    def __init__(self, args):
        self.ingtf = args.ingtf
        self.fasta = args.fasta
        self.species = args.species
        self.chrM = args.chrM
        self.genomeDir = args.genomeDir
        self.limitram = args.limitram
        self.threads = args.threads
        self.noindex = args.noindex

    def run(self):
        change_path()
        judgeFilexits(self.ingtf,self.fasta)
        if not self.noindex:
            star_index(self.fasta,self.ingtf,self.genomeDir,"%s/software/scStar"%__root_dir__,self.limitram,self.threads)
        write_config(self.genomeDir,self.species,self.chrM,self.fasta,self.ingtf)
        print("\033[0;32;40mAnalysis Complete\033[0m")

def mkref(args):
    Ref(args).run()

def helpInfo_mkref(parser):
    parser.add_argument(
        '--fasta',
        metavar='FASTA',
        help='Path to the FASTA file with the genome sequences.'
        )
    parser.add_argument(
        '--ingtf', 
        metavar='FILE' ,
        help='Path to the GTF file with annotations.'
        )
    parser.add_argument(
        '--species',
        metavar='STR',
        default='undefined',
        help='Species name, [default: undefined].'
        )
    parser.add_argument(
        '--chrM',
        metavar='STR',
        help='Mitochondrial chromosome name, [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--genomeDir',
        metavar='DIR',
        default=os.getcwd(),
        help='Path to the directory where genome files are stored, [default: current dir].'
        )
    parser.add_argument(
        '--limitram',
        metavar='INT',
        help='Maximum available RAM (bytes) for genome generation.',
        default=125000000000
        )
    parser.add_argument(
        '--threads',
        metavar='INT', 
        default=4,
        help='Number of threads used for the analysis, [default: 4].'
        )
    parser.add_argument(
        '--noindex',
        action='store_true',
        help='Only generate ref.json without constructing the database.'
        )
    return parser