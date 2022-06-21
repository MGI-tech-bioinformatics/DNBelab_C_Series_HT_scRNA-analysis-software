import os,glob,argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-o','--outdir', metavar='FILE', type=str,default=os.getcwd())
parser.add_argument('-n','--name', metavar='STR', type=str)
parser.add_argument('--combine', action='store_true', help="Combine sample metrics summary")
args = parser.parse_args()

    
def samplename():
    all_sample = []
    if args.name:
        for i in args.name.split(','):
            all_sample.append(i)
    else:
        html_list = glob.glob('%s/*/04.report/*.html'%(args.outdir))
        for html in html_list:
            sample = html.split('/')[-1].split('_')[0]
            all_sample.append(sample)
    return all_sample

def combine_result():
    result = os.path.join(args.outdir,'result')
    os.system('mkdir -p %s'%result)
    metrics_list = []
    html_list = []
    all_sample = samplename()
    for sample in all_sample:
        metrics_list.append('%s/%s/output/metrics_summary.xls'%(args.outdir,sample))
        html_list.append('%s/%s/output/%s_CDCPv2_scRNA_report.html'%(args.outdir,sample,sample))
        os.system('cp %s/%s/04.report/*.html %s'%(args.outdir,sample,result))
    metrics = pd.read_table("%s"%metrics_list[0],sep="\t",index_col=0)
    for i in range(1,len(metrics_list)):
        file_name = metrics_list[i]
        data = pd.read_table("%s"%file_name,sep="\t",index_col=0)
        metrics=pd.concat([metrics,data],axis=0)
    metrics.to_csv("%s/metrics_summary.xls"%result,sep='\t')

def remove_file():
    all_sample = samplename()
    os.system('rm -rf %s/Rplots.pdf'%(args.outdir))
    for sample in all_sample:
        os.system('rm -rf %s/%s/01.data/*_reads.fq.gz'%(args.outdir,sample))
        os.system('rm -rf %s/%s/01.data/*_barcode_counts_raw.txt'%(args.outdir,sample))
        os.system('rm -rf %s/%s/01.data/*.bam'%(args.outdir,sample))
        os.system('rm -rf %s/%s/02.count/*.bam'%(args.outdir,sample))
        os.system('rm -rf %s/%s/02.count/*_CB_UB_count.txt'%(args.outdir,sample))
        os.system('rm -rf %s/%s/02.count/cell_count_detail.xls'%(args.outdir,sample))

def main():
    if args.combine:
        combine_result()
        remove_file()
    else:
        remove_file()

if __name__=='__main__':
    main()
