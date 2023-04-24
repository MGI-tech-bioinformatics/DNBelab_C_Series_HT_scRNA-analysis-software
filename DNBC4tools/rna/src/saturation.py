#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import pysam
import multiprocessing
import time
import os
import shutil
import argparse
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline
from collections import defaultdict

import warnings
warnings.filterwarnings('ignore')

# Define command line arguments
parser = argparse.ArgumentParser(description='sequencing saturation')
parser.add_argument('-i','--inbam', metavar='FILE', type=str,help='input anno_decon_sorted.bam')
parser.add_argument('-f','--infile', metavar='FILE', type=str,help='input cellCount_report.csv')
parser.add_argument('-o','--outdir',help='storage outfile')
parser.add_argument('-q','--quality',type=int, default=20, help='Minimal map quality to filter. Default is 20')
parser.add_argument('-@','--threads',type=int, default=10, help='Analysis threads. Default is 10')
args = parser.parse_args()

# Extract values from command line arguments
inbam: str = args.inbam
infile: str = args.infile
outdir: str = args.outdir
quality: int = args.quality
threads: int = args.threads

def time_print(str):
    # Print time stamp along with message
    print("\033[32m%s\033[0m %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), str))

def count_result(in_bam:str, out_dir:str, ctg:str):
    # Count reads for a given contig (ctg) and output results to a file
    with open(os.path.join(out_dir,ctg+".txt"), 'w') as result:
        #result.write('\t'.join(['Cell', 'GeneID', 'UMI', 'Count']) + '\n')
        with pysam.AlignmentFile(in_bam, 'rb') as bam:
            gene_umi_dict = defaultdict(lambda:defaultdict(lambda: defaultdict(int)))
            for line in bam.fetch(contig=ctg):
                if int(line.flag) & 2304 != 0:
                    continue
                elif line.mapping_quality < quality or not line.has_tag('GN') :
                    continue
                elif len(line.get_tag('GN').split(';')) > 1:
                    continue
                else:
                    cell = line.get_tag('DB')
                    gene = line.get_tag('GN')
                    umi = line.get_tag('UB')
                gene_umi_dict[cell][gene][umi] += 1
            for cell in gene_umi_dict:
                for gene in gene_umi_dict[cell]:
                    for umi in gene_umi_dict[cell][gene]:
                        result.write('%s\t%s\t%s\t%s\n'%(cell,gene,umi,gene_umi_dict[cell][gene][umi]))

def combine_txt(indir:str):
    # Combine results from all contigs into a single file
    os.system("rm -rf %s/cell_count_detail.xls"%indir)
    #os.system("cat %s/temp/*.txt >> %s/cell_count_detail.xls"%(indir,indir))
    os.system("find %s/temp/ -name \"*.txt\" | xargs cat  >> %s/cell_count_detail.xls"%(indir,indir))

"""
This function calculates various saturation statistics for scRNA-seq data,
such as mean reads per cell, median genes per cell, sequencing saturation,
and UMI saturation, at different sampling fractions.
"""

def fraction_reads(outdir:str, totalreads:int, threads:int = 1):
    """
    Calculate saturation statistics for scRNA-seq data at different sampling fractions.

    Parameters:
        - outdir (str): Output directory for the results.
        - totalreads (int): Total number of reads in the scRNA-seq dataset.
        - threads (int): Number of threads to use for parallel processing (default: 1).

    Returns:
        None. The results are saved in a tab-separated file called "saturation.xls" in the
        specified output directory.
    """

    # Read the cell count detail file and convert it to a Polars dataframe.
    cellcount_df = pl.read_csv(
        os.path.join(outdir, "cell_count_detail.xls"),
        has_headers=False,
        sep='\t',
        use_pyarrow=False,
        n_threads=threads,
        columns=["column_1", "column_2", "column_3", "column_4"],
        new_columns=["Cell", "GeneID", "UMI", "Count"]
    ).with_columns([
        pl.col("Cell").cast(pl.Categorical),
        pl.col("GeneID"),
        pl.col("UMI").cast(pl.Categorical),
        pl.col("Count").cast(pl.UInt32)
    ])

    # Define the sampling fractions to use.
    sampling_fractions = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    sampling_fractions_length = len(sampling_fractions)

    # Create an empty Pandas dataframe to store the results.
    stats_df = pd.DataFrame(
        {
            "Mean Reads per Cell": np.zeros(sampling_fractions_length, np.uint32),
            "Median Genes per Cell": np.zeros(sampling_fractions_length, np.uint32),
            "Mean Gene per Cell": np.zeros(sampling_fractions_length, np.uint32),
            "Total Gene": np.zeros(sampling_fractions_length, np.uint32),
            "Sequencing Saturation": np.zeros(sampling_fractions_length, np.float64),
            "UMI Saturation": np.zeros(sampling_fractions_length, np.float64)
        },
        index=pd.Index(data=np.array(sampling_fractions), name="sampling_fraction")
    )

    # Calculate the saturation statistics for the full dataset.
    cellcount_all_df = cellcount_df.with_column(pl.col("Count").repeat_by(pl.col("Count"))).explode("Count")
    stats_df.loc[1.0, "Mean Reads per Cell"] = round(totalreads/(pl.n_unique(cellcount_df['Cell'])))
    stats_df.loc[1.0, "UMI Saturation"] = round((1- (cellcount_df.filter(pl.col("Count") == 1).height)/(cellcount_df.height))*100,2)
    stats_df.loc[1.0, "Sequencing Saturation"] = round((1- cellcount_df.height/cellcount_df.filter(pl.col('Cell')!='None').select([pl.col("Count").sum()])[0,0])*100,2)
    cellcount_df = cellcount_df.with_column(pl.col("GeneID").str.split(";")).explode("GeneID")
    stats_df.loc[1.0, "Median Genes per Cell"] = round(cellcount_df.filter(pl.col('Cell')!='None').groupby("Cell").agg(
            [pl.n_unique(['GeneID']).alias("MedianGene")]).select([pl.col("MedianGene").median()])["MedianGene"][0])
    stats_df.loc[1.0, "Mean Gene per Cell"] = round(cellcount_df.filter(pl.col('Cell')!='None').groupby("Cell").agg(
            [pl.n_unique(['GeneID']).alias("MeanGene")]).select([pl.col("MeanGene").mean()])["MeanGene"][0])
    stats_df.loc[1.0, "Total Gene"] = pl.n_unique(cellcount_df['GeneID'])
    del cellcount_df

    for sampling_fraction in sampling_fractions:
        if sampling_fraction == 0.0:
            continue
        elif sampling_fraction == 1.0:
            continue
        else:
            cellcount_sampled=cellcount_all_df.sample(frac=sampling_fraction)
            cellcount_sampled=cellcount_sampled.groupby(["Cell", "GeneID","UMI"]).agg([pl.col("UMI").count().alias("Count")])
            stats_df.loc[sampling_fraction, "Mean Reads per Cell"] = round(totalreads*float(sampling_fraction)/pl.n_unique(cellcount_sampled['Cell']))
            stats_df.loc[sampling_fraction, "UMI Saturation"] = round((1- (cellcount_sampled.filter(pl.col("Count") == 1).height)/cellcount_sampled.height)*100,2)
            stats_df.loc[sampling_fraction, "Sequencing Saturation"] = round((1- (cellcount_sampled.height-1)/cellcount_sampled.filter(pl.col('Cell')!='None').select([pl.col("Count").sum()])[0,0])*100,2)
            cellcount_sampled = cellcount_sampled.with_column(pl.col("GeneID").str.split(";")).explode("GeneID")
            stats_df.loc[sampling_fraction, "Median Genes per Cell"] = round(cellcount_sampled.filter(pl.col('Cell')!='None').groupby("Cell").agg(
                    [pl.n_unique(['GeneID']).alias("MedianGene")]).select([pl.col("MedianGene").median()])["MedianGene"][0])
            stats_df.loc[sampling_fraction, "Mean Gene per Cell"] = round(cellcount_sampled.filter(pl.col('Cell')!='None').groupby("Cell").agg(
                    [pl.n_unique(['GeneID']).alias("MeanGene")]).select([pl.col("MeanGene").mean()])["MeanGene"][0])
            stats_df.loc[sampling_fraction, "Total Gene"] = pl.n_unique(cellcount_sampled['GeneID'])
            del cellcount_sampled
    del cellcount_all_df
    stats_df.to_csv(os.path.join(outdir,"saturation.xls"),sep='\t')

def main():
    df = pd.read_csv(open(infile),encoding="utf_8",dtype=str,header=None,sep=",")
    totalReads = int(df[1][2])
    resulttemp= os.path.join(outdir,"temp")
    if not os.path.exists(resulttemp):
        os.system('mkdir -p %s'%resulttemp)
    time_print("count reads")
    with pysam.AlignmentFile(inbam, 'rb') as bam:
        ctg_list = []
        for item in bam.header["SQ"]:
            item = dict(item)
            ctg_list.append(item['SN'])
    pool1 = multiprocessing.Pool(processes=threads)
    for ctg in ctg_list:
        pool1.apply_async(count_result,(inbam, resulttemp, ctg,))
    pool1.close()
    pool1.join()
    time_print("combine txt")
    combine_txt(outdir)
    shutil.rmtree(resulttemp)
    time_print("Analysis saturation")
    fraction_reads(outdir,totalReads)
    time_print("plot saturation")
    plot_saturation()
    time_print("Analysis complete")

def to_percent(temp: float, position: int) -> str:
    """
    Convert a number to a percentage string with 'k' added to the end if over 1000.

    Args:
    temp: The number to be converted.
    position: Not used.

    Returns:
    The percentage string with 'k' added if necessary.
    """
    return '%d'%(temp/1000) + 'k'

def umi_saturation(ax: plt.axes, table: pd.DataFrame):
    """
    Plot the UMI sequencing saturation curve.

    Args:
    ax: The plot axes to draw the curve on.
    table: The data table with columns 'Mean Reads per Cell' and 'Sequencing Saturation'.
    """
    xnew = np.linspace(table['Mean Reads per Cell'].min(),table['Mean Reads per Cell'].max(),300)
    smooth = make_interp_spline(table['Mean Reads per Cell'],table['Sequencing Saturation']/100)(xnew)
    ax.set_xlim([0, table['Mean Reads per Cell'].max()])
    ax.set_ylim([0, 0.9999])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(to_percent))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False)
    ax.grid(linestyle='--')
    ax.plot(xnew,smooth,linewidth=3.0)
    ax.axhline(y=0.9,ls="--",c="black",linewidth=2.0)
    ax.set(xlabel='Mean Reads per Cell', ylabel='Sequencing Saturation',title='Sequencing Saturation')

def gene_saturation(ax: plt.axes, table: pd.DataFrame):
    """
    Plot the gene expression saturation curve.

    Args:
    ax: The plot axes to draw the curve on.
    table: The data table with columns 'Mean Reads per Cell' and 'Median Genes per Cell'.
    """
    xnew = np.linspace(table['Mean Reads per Cell'].min(),table['Mean Reads per Cell'].max(),300)
    smooth = make_interp_spline(table['Mean Reads per Cell'],table['Median Genes per Cell'])(xnew)
    ax.set_xlim([0, table['Mean Reads per Cell'].max()])
    ax.set_ylim([0, table['Median Genes per Cell'].max()])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(to_percent))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False)
    ax.grid(linestyle='--')
    ax.plot(xnew,smooth,linewidth=3.0)
    ax.set(xlabel='Mean Reads per Cell', ylabel='Median Gene per Cell',title='Median Gene per Cell')

def plot_saturation():
    """
    Generate and save the UMI and gene expression saturation plots.
    """
    for_plot = pd.read_table(os.path.join(outdir,'saturation.xls'),sep='\t')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), tight_layout=True)
    arts = umi_saturation(ax1,for_plot)
    arts = gene_saturation(ax2,for_plot)
    fig.savefig(os.path.join(outdir,'saturation.png'),facecolor='white',transparent=False,dpi=400)
    plt.close(fig)


if __name__=='__main__':
    main()
