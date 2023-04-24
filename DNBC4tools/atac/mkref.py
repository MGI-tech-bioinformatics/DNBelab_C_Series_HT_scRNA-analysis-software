import os, re, json, collections
from dnbc4tools.tools.utils import judgeFilexits, change_path, bin_path, start_print_cmd
from dnbc4tools.__init__ import __root_dir__
from dnbc4tools.tools.mkgtf import read_gtf

# Function to create a reference configuration file for genome processing
def write_config(genomeDir, species, gtf, chrM, blacklist, fasta,tag,prefix):

    # Fetching path for samtools
    bin_command = bin_path()

    # Creating directory for storing genome index files
    if not os.path.exists(genomeDir):
        os.system('mkdir -p %s' % genomeDir)

    # Creating fasta index file using samtools
    faidx_cmd = '%s/samtools faidx -o %s/genome.fa.fai %s' % (bin_command, genomeDir, fasta)
    start_print_cmd(faidx_cmd)

    # Getting genome size using helper functions
    genome_size = get_genome_size(genomeDir, species, gtf,prefix,tag)

    # Creating a dictionary to store reference file locations and information
    ref_dict = collections.OrderedDict()
    ref_dict['species'] = str(species)
    ref_dict['genome'] = os.path.abspath(fasta)
    ref_dict['index'] = os.path.abspath('%s/genome.index' % genomeDir)
    ref_dict['chromeSize'] = os.path.abspath('%s/chrom.sizes' % genomeDir)
    ref_dict['tss'] = os.path.abspath('%s/tss.bed' % genomeDir)
    ref_dict['promoter'] = os.path.abspath('%s/promoter.bed' % genomeDir)

    # Checking if blacklist file is provided
    if blacklist != 'None':
        if blacklist == 'hg19':
            ref_dict['blacklist'] = '%s/config/blacklist/hg19.full.blacklist.bed'%__root_dir__
        elif blacklist == 'hg38':
            ref_dict['blacklist'] = '%s/config/blacklist/hg38.full.blacklist.bed'%__root_dir__
        elif blacklist == 'mm10':
            ref_dict['blacklist'] = '%s/config/blacklist/mm10.full.blacklist.bed'%__root_dir__
        else:
            if os.path.exists(blacklist):
                ref_dict['blacklist'] = os.path.abspath(blacklist)
            else:
                raise Exception('Unable to find blacklist file!')
    else:
        ref_dict['blacklist'] = 'None'

    # Getting mitochondrial chromosome information using helper function
    ref_dict['chrmt'] = get_chrM(chrM,genomeDir)

    # Adding genome size information to the dictionary
    ref_dict['genomesize'] = str(genome_size)

    # Writing reference dictionary to a JSON file
    with open(('%s/ref.json' % genomeDir), 'w', encoding='utf-8') as (jsonfile):
        json.dump(ref_dict, jsonfile, indent=4, ensure_ascii=False)
        jsonfile.write('\n')

# Function to get mitochondrial chromosome information
def get_chrM(chrM,genomeDir):
    """
    Given a mitochondrial chromosome name or 'auto', return the correct name of the mitochondrial chromosome in a genome.

    :param chrM (str): Name of the mitochondrial chromosome or 'auto' to automatically detect the name.
    :param genomeDir (str): Path to directory containing genome files.
    :return chrMT (str): Correct name of the mitochondrial chromosome in the genome.
    """
    if chrM == 'auto':

        # List of possible names for mitochondrial chromosome
        chrmtlist = ['chrM','MT','chrMT','mt','Mt']
        
        # Checking if fasta index file exists
        if os.path.exists(os.path.join(genomeDir,'genome.fa.fai')):
            with open(os.path.join(genomeDir,'genome.fa.fai'),'r') as file:
                
                # Extracting chromosome names from fasta index file
                chrlist =  [line.split()[0] for line in file]
                
                # Finding intersection between possible names and actual chromosome names
                union = list(set(chrmtlist) & set(chrlist))
                if union:
                    chrMT = union[0]
                else:
                    chrMT = 'None'
        else:
            chrMT = 'None'
    else:
        chrMT = chrM
    print('Mitochondrial chromosome: %s'%chrMT)
    return chrMT

def get_genome_size(genomeDir,species,gtf,prefix, tag):
    """
    Calculate the total size of the genome and return a species-specific identifier.

    Parameters:
    - genomeDir: directory path where genome files are stored.
    - species: name of the species (human, mouse, etc.)
    - gtf: file path to the GTF annotation file.
    - prefix: optional prefix for chromosome names.
    - tag: tag used to label TSS promoter regions.

    Returns:
    - genome_size: string identifier for the genome size.

    """

    # Open the fasta index file and write chromosome sizes to a file
    with open('%s/genome.fa.fai'%genomeDir,'r') as fai:
        with open('%s/chrom.sizes'%genomeDir,'w') as result:
            
            # Get a list of primary contigs from the GTF annotation
            primary_contigs = get_primary_contigs(gtf)
            
            # Filter chromosomes by prefix and primary contigs
            filtered_chromosomes = filter_chromosomes(primary_contigs, prefix)
            
            # Print out the remaining chromosomes
            print('The remaining chromosomes are:')
            print(filtered_chromosomes)
            
            # Label TSS promoter regions for each gene
            get_tss_promoter(gtf, tag, genomeDir,filtered_chromosomes)
            
            # Iterate over the fasta index file and extract chromosome sizes
            chr_size = []
            for line in fai:
                line = line.strip()
                lst = line.split('\t')
                
                # Only include filtered chromosomes
                if lst[0] in filtered_chromosomes:
                    result.write(lst[0]+'\t'+lst[1]+'\n')
                    chr_size.append(int(lst[1]))
    
    # Calculate the total genome size
    genome_size = sum(chr_size)

    # Assign a species-specific identifier based on the species name
    if species in  ['Human','hg19','hg38', 'Homo_sapiens','GRCh38', 'GRCh37']:
        genome_size = 'hs'
    elif species in ['Mouse', 'mm10' ,'Mus_musculus','GRChm38']:
        genome_size = 'mm'
    elif species == 'Caenorhabditis elegans':
        genome_size = 'ce'
    elif species == 'Fruitfly':
        genome_size = 'dm'
    else:
        genome_size = genome_size
    return genome_size

def get_tss_promoter(gtf, type, genomeDir,chromosome):
    """
    This function extracts the TSS (Transcription Start Site) and promoter regions from a GTF file.
    The TSS is defined as the first base of the transcript (either start or end based on strand).
    The promoter region is defined as 2kb upstream and 2kb downstream of the TSS.

    Parameters:
    - genomeDir - directory to write the output BED files
    - chromosome - a list of chromosomes to filter by (optional)
    - gtf - path to the GTF file
    - type - type of annotation to extract (e.g. 'exon', 'gene', 'transcript')

    """

    # Read the GTF file and filter by the given type
    gtfread = read_gtf(gtf, type)
    
    # Initialize empty sets to store TSS and promoter regions
    tss_bed = set()
    promoter_bed = set()
    
    # Initialize a flag for whether any transcript has the 'basic' tag
    no_basic = True
    
    # Open files to write the output BED files for TSS and promoter regions
    tssfile = open(('%s/tss.bed' % genomeDir), 'w', encoding='utf-8')
    promoterfile = open(('%s/promoter.bed' % genomeDir), 'w', encoding='utf-8')

    # Loop through each line in the filtered GTF file
    for type_cell in gtfread:
        if type_cell:
            if type_cell[0].startswith('#'):
                pass
            else:
                cell_lst = type_cell[0].split('\t')
                chrome = cell_lst[0]
                if chromosome and chrome not in chromosome:
                    continue

                # Determine the TSS and promoter regions based on strand
                start = int(cell_lst[3])
                end = int(cell_lst[4])
                strand = cell_lst[6]
                if strand == '+':
                    tss_start = start - 1
                    tss_end = tss_start + 1
                    promoter_start = start - 1 - 2000
                    promoter_end = start + 2000
                elif strand == '-':
                    tss_start = end - 1
                    tss_end = tss_start + 1
                    promoter_start = end - 1 - 2000
                    promoter_end = end + 2000

                # Parse the transcript tags and extract the gene name or ID
                transcript_tags = re.findall('tag\\s["*](\\w+)*["*]', cell_lst[-1])
                pattern = re.compile('(\\S+?)\\s*"(.*?)"')
                aDict = collections.OrderedDict()
                for m in re.finditer(pattern, cell_lst[-1]):
                    key = m.group(1)
                    value = m.group(2)
                    aDict[key] = value
                else:
                    if 'gene_name' in aDict:
                        gene_id = aDict['gene_name']
                    else:
                        if 'gene_id' in aDict:
                            gene_id = aDict['gene_id']
                        else:
                            gene_id = ''
                    is_basic = 'basic' in transcript_tags
                    if is_basic:
                        no_basic = False

                    # Create tuples for the TSS and promoter regions and add them to the respective sets
                    tss = (str(chrome), str(tss_start), str(tss_end), str(gene_id), '.', str(strand), is_basic)
                    promoter = (str(chrome), str(promoter_start), str(promoter_end), str(gene_id), is_basic)
                    tss_bed.add(tss)
                    promoter_bed.add(promoter)

    for tss in sorted(tss_bed, key=(lambda x: (x[0], int(x[1])))):
        if no_basic or tss[-1]:
            tssfile.write('\t'.join(tss[0:-1]) + '\n')
    for promoter in sorted(promoter_bed, key=(lambda x: (x[0], int(x[1])))):
        if no_basic or promoter[-1]:
            promoterfile.write('\t'.join(promoter[0:-1]) + '\n')
    tssfile.close()
    promoterfile.close()


def get_primary_contigs(gtf):
    """
    This function extracts the primary contigs from a GTF file.

    Args:
    - gtf (str): The path to the GTF file.

    Returns:
    - primary_contigs (set): A set of primary contigs extracted from the GTF file.
    """
    contigs = []
    with open(gtf, 'r') as (gtffile):
        for line in gtffile:
            line = line.strip()
            if line.startswith('#'):
                pass
            elif line == '':
                pass
            else:
                lst = line.split('\t')
                if lst[2] == 'gene' or lst[2] == 'transcript':
                    contigs.append(lst[0])
    primary_contigs = set(contigs)
    return primary_contigs

def filter_chromosomes(chromosomes, prefix=None):
    """
    Filter chromosomes by prefix or full name.

    :param chromosomes: A list of chromosome names.
    :param prefix: A string, a list of strings, or a comma-separated string representing the prefix(es) or full name(s) of the chromosome(s) to keep.
    :return: A list of chromosome names that match the given prefix(es) or full name(s).
    """
    if prefix is None:
        return chromosomes
    elif isinstance(prefix, str):
        prefix = [p.strip() for p in prefix.split(",")]
    elif isinstance(prefix, list):
        prefix = [p.strip() for p in prefix]
    else:
        raise TypeError("prefix must be a string, a list of strings, or a comma-separated string")

    filtered_chromosomes = []
    for chromosome in chromosomes:
        if any(chromosome.startswith(p) or chromosome == p for p in prefix):
            filtered_chromosomes.append(chromosome)
    
    if not filtered_chromosomes:
        return chromosomes

    return filtered_chromosomes


def atac_index(fasta, genomeDir):
    index_cmd = '%s/software/chromap --build-index --ref %s --output %s/genome.index' % (__root_dir__, fasta, genomeDir)
    print('chromap verison: 0.2.3-r407')
    print('runMode: genomeGenerate')
    print('genomeDir: %s' % os.path.abspath(genomeDir))
    start_print_cmd(index_cmd)


class Ref:
    """
    This class represents the reference genome and annotation information required for ATAC-seq data analysis.

    Args:
        args (argparse.Namespace): Namespace containing the command-line arguments.

    Attributes:
        ingtf (str): Path to the GTF file with annotations.
        fasta (str): Path to the FASTA file with the genome sequences.
        genomeDir (str): Path to the directory where genome files are stored.
        species (str): Species name.
        tag (str): Select the type to generate BED file.
        chrM (str): Mitochondrial chromosome name.
        blacklist (str): Path to a file containing genomic regions that exhibit high signal noise.
        prefix (str): Filter chromosomes by prefix or full name.
        noindex (bool): Whether or not to construct the database.

    Methods:
        run(): Runs the reference genome and annotation information generation process.
    """

    def __init__(self, args):
        self.ingtf = args.ingtf
        self.fasta = args.fasta
        self.genomeDir = args.genomeDir
        self.species = args.species
        self.tag = args.tag
        self.chrM = args.chrM
        self.blacklist = args.blacklist
        self.prefix = args.prefix
        self.noindex = args.noindex

    def run(self):
        """
        Runs the reference genome and annotation information generation process. 
        
        Changes the working directory, checks if the input files exist, builds an index if necessary, 
        writes configuration files, and prints a completion message.
        """

        change_path()
        judgeFilexits(self.ingtf, self.fasta)
        print('\x1b[0;32;40mBuilding index for dnbc4tools atac\x1b[0m')
        if not self.noindex:
            atac_index(self.fasta, self.genomeDir)
        write_config(self.genomeDir, self.species, self.ingtf, self.chrM, self.blacklist, self.fasta, self.tag, self.prefix)
        print('\x1b[0;32;40mAnalysis Complete\x1b[0m')


def mkref(args):
    """
    A function that creates a Ref object and runs the reference genome and annotation information generation process.
    
    Args:
        args (argparse.Namespace): Namespace containing the command-line arguments.
    """

    Ref(args).run()


def helpInfo_mkref(parser):
    """
    Adds arguments to an argparse parser for creating the reference genome and annotation information.
    
    Args:
        parser (argparse.ArgumentParser): An ArgumentParser object.
        
    Returns:
        argparse.ArgumentParser: The modified ArgumentParser object.
    """
    
    parser.add_argument(
        '--fasta',
        metavar='FASTA',
        help='Path to the fasta file with the genome sequences.'
        )
    parser.add_argument(
        '--ingtf',
        metavar='FILE',
        help='Path to the gtf file with annotations.'
        )
    parser.add_argument(
        '--genomeDir',
        metavar='DIR',
        default=(os.getcwd()),
        help='Path to the directory where genome files are stored, [default: current dir].'
        )
    parser.add_argument(
        '--species',
        metavar='STR',
        default='undefined',
        help='Species name, [default: undefined].'
        )
    parser.add_argument(
        '--tag',
        metavar='STR',
        default='transcript',
        help='Select the type to generate bed, [default: transcript].'
        )
    parser.add_argument(
        '--chrM',
        metavar='STR',
        default='auto',
        help='Mitochondrial chromosome name, [default: auto].'
        )
    parser.add_argument(
        '--blacklist',
        metavar='FILE',
        default='None',
        help='Genomic regions that are known to exhibit high signal noise, [default: None].'
        )
    parser.add_argument(
        '--prefix',
        default='None',
        help='Filter chromosomes by prefix or full name, [default: None].'
        )
    parser.add_argument(
        '--noindex',
        action='store_true',
        help='Only generate ref.json without constructing database.'
        )
    return parser