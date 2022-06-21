import argparse,os
parser = argparse.ArgumentParser(description='Create star index')
parser.add_argument('--action', metavar='SELECT', default='mkref', choices=['mkref', 'mkgtf', 'stat'], help='Select the action for your program, include mkref,mkgtf,stat, default is mkref') 
parser.add_argument('--ingtf', metavar='FILE' ,help='Set ingtf in mkref,mkgtf or stat')
parser.add_argument('--outgtf',metavar='FILE', help='Set outgtf in mkgtf')
parser.add_argument('--attribute',metavar='DICT',default = ['gene_type=protein_coding'], nargs='+',help='Set the filter parameter in mkgtf, Key-value pair in attributes field to be kept in the GTF, \
    default is gene_type:protein_coding, you can set up multiple groups separated by blankspace')
parser.add_argument('--outstat',metavar='FILE', default = 'gtf_type.txt',help='Set the stats outfile in stat, default is "gtf_type.txt" in current dir')
parser.add_argument('--type',metavar='STR', default = 'gene_type',help='Set the the type for stat, default is gene_type')
parser.add_argument('--fasta',metavar='FASTA',help='Set the fasta in mkref')
parser.add_argument('--star_dir',metavar='DIR',default=os.getcwd(), help='Set the star indexdir in mkref, default is current dir')
parser.add_argument('--star',metavar='PROGRAM',help='Set the star program path in mkref')
parser.add_argument('--thread',metavar='INT', help='Set the threads in mkref')
args = parser.parse_args()

from typing import List, Iterable
from collections import Counter
from subprocess import check_call

def read(fp: str, feature: str) -> Iterable[List[str]]:
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

def filtergtf(gtf,filtergtf,attribute):
    d = dict(i.split(":")[::-1] for i in attribute)
    gtfread = read(gtf,'gene')
    result = open(filtergtf,'w')
    for i in gtfread:
        if i[0].startswith('#'):
            result.writelines(i)
        else:
            lst = i[0].split('\t')[-1].replace('"', '').strip().strip(';')
            lstlist = lst.split('; ')
            aDict = {}
            for j in lstlist:
                aDict[j.split(' ')[0]] = j.split(' ')[-1]
            for key1,value1 in aDict.items():
                for key2,value2 in d.items():
                    if key1 == value2 and key2 == value1:
                        result.writelines(i)
    result.close()

def statgtf(gtf,keyword,outfile):
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
                    lstlist = lst[-1].replace('"', '').strip().strip(';')
                    llstlist = lstlist.split('; ')
                    aDict = dict(l.split(' ') for l in llstlist)
                    sumDict.append(aDict[keyword])
        result = Counter(sumDict)
        outfile = open(outfile,'w')
        outfile.write('Type'+'\t'+'Count'+'\n')
        for k,v in sorted(result.items(), key = lambda x:x[1], reverse=True):
            outfile.write(f'{k}\t{v}\n')
        outfile.close()

def star_index(fasta,gtf,star_dir,star_program,threads):
    if not os.path.exists(star_dir):
        os.system('mkdir -p %s'%star_dir)
    star_cmd = '%s --runMode genomeGenerate --runThreadN %s --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang 99 --limitGenomeGenerateRAM 125000000000'\
        %(star_program,threads,star_dir,fasta,gtf)
    print('STAR verison: 2.7.2b')
    print('runMode: genomeGenerate')
    print('runThreadN: %s'%threads)
    print('genomeDir: %s'%star_dir)
    print('fasta: %s'%fasta)
    print('gtf: %s'%gtf)
    check_call(star_cmd,shell=True)

if __name__=='__main__':
    if args.action == 'stat':
        statgtf(args.ingtf,args.type,args.outstat)
    if args.action == 'mkgtf':
        filtergtf(args.ingtf,args.outgtf,args.attribute)
    if args.action == 'mkref':
        star_index(args.fasta,args.ingtf,args.star_dir,args.star,args.thread)

