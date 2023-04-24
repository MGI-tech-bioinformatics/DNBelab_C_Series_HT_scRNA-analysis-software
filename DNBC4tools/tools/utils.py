import os,sys,json
import time
import logging
import sys
from datetime import timedelta
from subprocess import check_call
from dnbc4tools.__init__ import __root_dir__

def str_mkdir(arg):
    if not os.path.exists(arg):
        os.system('mkdir -p %s'%arg)

def change_path():
    os.environ['PATH'] += ':'+'/'.join(str(__root_dir__).split('/')[0:-4])+ '/bin'
    os.environ['LD_LIBRARY_PATH'] = '/'.join(str(__root_dir__).split('/')[0:-4]) + '/lib'

def bin_path():
    bin_command = '/'.join(str(__root_dir__).split('/')[0:-4])+ '/bin'
    return bin_command
    
def rm_temp(*args):
    for filename in args:
        if os.path.exists(filename):
            os.remove(filename)
        else:
            pass

def start_print_cmd(arg):
    #print(arg)
    check_call(arg,shell=True)

def logging_call(popenargs,name,dir):
    today = time.strftime('%Y%m%d', time.localtime(time.time()))
    logfile = '%s/log/%s.txt'%(dir,today)
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(level = logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.ERROR)
        console_handler.formatter = formatter
        file_handler = logging.FileHandler(logfile,encoding="utf8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

    logger.info('Promgram start...')
    logger.info(popenargs)
    start = time.time()
    check_call(popenargs,shell=True)
    #logger.info(output.decode('utf-8'))
    logger.info('Promgram end...')
    end = time.time()
    used = timedelta(seconds=end - start)
    logger.info('Program time used: %s', used)
    logger.info('\n')  

def judgeFilexits(*args):
    for input_files in args:
        for input_file in input_files.split(','):
            if not os.path.exists(input_file): 
                print(" ------------------------------------------------") 
                print("Error: Cannot find input file or dir %s."%(str(input_file))) 
                print(" ------------------------------------------------") 
                sys.exit()
            else:
                pass

def hamming_distance(chain1, chain2):
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(chain1, chain2))))


def read_json(file):
    with open(file,'r',encoding='utf8')as fp:
        json_data = json.load(fp)
    return json_data

def seq_comp(seq):
    nt_comp = {'A':'0', 'C':'1', 'G':'2', 'T':'3'}
    length = len(seq)-1
    sum = 0
    for k,v in enumerate(seq.upper()):
        sum += int(nt_comp[v])*(4**(length-k))
    return str('%010x'%sum).upper()

def read_anndata(path):
    import scipy.io
    from scipy.sparse import csr_matrix
    import anndata
    import pandas as pd
    mat = scipy.io.mmread(path+"/"+"matrix.mtx.gz").astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat,dtype="float32")
    genes = pd.read_csv(path+'/'+'features.tsv.gz', header=None, sep='\t')
    var_names = genes[0].values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    adata.obs_names = pd.read_csv(path+'/'+'barcodes.tsv.gz', header=None)[0].values
    adata.var_names_make_unique()
    return adata

def png_to_base64(file,filename,outdir):
    import base64
    base64_path = outdir +'/'+filename+'.base64'
    if os.path.isfile(file):
        with open(file, "rb") as f:
            base64_data = base64.b64encode(f.read())
            s = base64_data.decode()
            base64_path_f = open(base64_path, 'w')
            base64_path_f.write('<img src=data:image/'+'png'+';base64,'+s+">")
            base64_path_f.close()

def csv_datatable(file,outfile):
    import pandas as pd
    if os.path.exists(file):
        df= pd.read_csv(open(file),encoding="utf-8",dtype=str,)
        fw = open(outfile,'w')
        for index, row in df.iterrows():
            fw.write('<tr><td>'+row['gene']+'</td>'\
                +'<td>'+row['cluster']+'</td>'\
                +'<td>'+row['p_val_adj']+'</td>'\
                +'<td>'+row['p_val']+'</td>'\
                +'<td>'+row['avg_log2FC']+'</td>'\
                +'<td>'+row['pct.1']+'</td>'\
                +'<td>'+row['pct.2']+'</td>'\
            )
        fw.close()

def write_matrix(adata,outdir):
    import scipy.io
    import scipy.sparse
    import shutil
    import gzip
    adata.X = scipy.sparse.csr_matrix(adata.X.astype('int32'))
    scipy.io.mmwrite(
        "%s/matrix.mtx"%outdir, 
        adata.X.transpose()
        )
    adata.var.to_csv(
        '%s/features.tsv.gz'%outdir, 
        sep='\t', index=True, header=False,
        compression='gzip'
        )
    adata.obs.to_csv(
        '%s/barcodes.tsv.gz'%outdir, 
        sep='\t', index=True, header=False,
        compression='gzip'
        )
    with open("%s/matrix.mtx"%outdir,'rb') as mtx_in:
        with gzip.open("%s/matrix.mtx"%outdir + '.gz','wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove("%s/matrix.mtx"%outdir)