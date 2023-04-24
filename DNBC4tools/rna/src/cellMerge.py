#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import argparse
import collections
import pandas as pd
from itertools import groupby
from typing import List, Dict
from dnbc4tools.tools.utils import seq_comp

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def barcodeTranslatefile(combined_file: str, select_barcode: str, barcodeTranslate: str,
                         barcodeTranslate_hex: str, cellid: str) -> None:
    """
    The function reads in the combined_file and extract the barcode and cell names of beads
    with low cell number, storing them in a dictionary. It then reads in the select_barcode file,
    assigning new cell names to barcodes that are not in the dictionary. The resulting dictionary
    is then written to barcodeTranslate and barcodeTranslate_hex files, and the unique cell names
    are written to cellid file.

    Args:
    - combined_file: str, path of the input combined file
    - select_barcode: str, path of the input select barcode file
    - barcodeTranslate: str, path of the output barcode translate file
    - barcodeTranslate_hex: str, path of the output barcode translate hex file
    - cellid: str, path of the output cell id file

    Returns:
    None
    """
    multi_barcodelist: List[str] = []
    multi_celllist: List[str] = []
    multibeads_filter_dict: dict = {}
    one_barcodedict: dict = collections.OrderedDict()

    with open(combined_file,'r') as multi_beads:
        for line in multi_beads:
            line = line.strip()
            if line:
                multi_barcodelist.append(line.split('\t')[0])
                multi_celllist.append(line.split('\t')[-1])
                if int(line.split('\t')[-1].split('_N')[-1]) <= 9:
                    multibeads_filter_dict.setdefault(line.split('\t')[-1],[]).append(line.split('\t')[0])

    with open(select_barcode,'r') as allbarcode:
        if multi_celllist:
            num = int(multi_celllist[-1].split('CELL')[-1].split('_')[0])
        else:
            num = 0
        for line in allbarcode:
            line = line.strip()
            if line not in multi_barcodelist:
                num += 1
                cellName = 'CELL'+str(num)+'_N1'
                one_barcodedict[line] = cellName

    with open(barcodeTranslate,'w') as barcode_cell,\
        open(barcodeTranslate_hex,'w') as barcode_cell_hex,\
        open(cellid,'w') as cellFile:
        celllist = []
        for k2,v2 in one_barcodedict.items():
            barcode_cell.write(f'{k2}\t{v2}'+'\n')
            celllist.append(v2)
            k2_hex = seq_comp(k2)
            barcode_cell_hex.write(f'{k2_hex}\t{v2}'+'\n')
        for k1,v1 in multibeads_filter_dict.items():
            for v1terms in v1:
                barcode_cell.write(f'{v1terms}\t{k1}'+'\n')
                celllist.append(k1)
                v1terms_hex = seq_comp(v1terms)
                barcode_cell_hex.write(f'{v1terms_hex}\t{k1}'+'\n')
        cellRmDup = [x for x, _ in groupby(celllist)]
        cellFile.write('\n'.join(cellRmDup))

def merge_graph(figtable: pd.DataFrame, cellnum: int, outdir: str) -> None:
    """
    Merge multiple bar graphs into one and save as png and pdf files.

    Args:
    - figtable: A pandas DataFrame containing the data to plot.
    - cellnum: An integer representing the total number of cells.
    - outdir: A string specifying the output directory.

    Returns:
    - None
    """
    # Set plot parameters
    params: Dict[str, object] = {
        'figure.figsize': (7.65, 5.72),
        'axes.labelsize': 'larger',
        'figure.dpi': 100,
        'axes.titlelocation': 'left',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.handlelength': 1.2,
        'legend.handleheight': 1.2,
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium'
    }
    plt.rcParams.update(params)

    # Set color palette
    set2_colors: List[str] = sns.color_palette('Set2', n_colors=len(figtable))

    # Plot bar graph
    ax = sns.barplot(x='Num', y='Count', data=figtable, palette=set2_colors, saturation=1)

    # Add legend
    legend_labels: List[str] = []
    for i in range(len(figtable)):
        label = figtable.iloc[i]['Num'] + " " + str(figtable.iloc[i]['Count'])
        legend_labels.append(label)

    handles: List[Patch] = [
        Patch(facecolor=set2_colors[i], edgecolor='none', label=legend_labels[i])
        for i in range(len(figtable))
    ]
    ax.legend(handles=handles, loc='upper right', frameon=False, labelspacing=1)

    # Set axis labels and title
    ax.set(xlabel='Number of beads per droplet', ylabel='CellCount',
           title='Total cell number %s' % cellnum)

    # Save plot as png and pdf files
    plt.savefig('%s/cellNumber_merge.png' % outdir)
    plt.savefig('%s/cellNumber_merge.pdf' % outdir)

def summary_count(combined_file: str, select_barcode: str, beads_stat: str, barcodeTranslate: str,
                  barcodeTranslate_hex: str, cellid: str, cellCount: str, outdir: str) -> None:
    """
    Generate a count report of cell and barcode summary statistics.
    This function takes several file paths and an output directory as input, 
    and produces a summary report and a figure as output.
    Call the barcodeTranslatefile function to translate the barcodes in the combined file 
    and save the output to the barcodeTranslate and barcodeTranslate_hex files.

    Args:
        combined_file (str): Path to the combined list file for analysis.
        select_barcode (str): Path to the selected barcode file.
        beads_stat (str): Path to the beads stat file.
        barcodeTranslate (str): Path to the barcode translate file.
        barcodeTranslate_hex (str): Path to the barcode translate hex file.
        cellid (str): Path to the cell id file.
        cellCount (str): Path to the cell count report file.
        outdir (str): Path to the output directory for analysis.

    """
    barcodeTranslatefile(combined_file, select_barcode, barcodeTranslate, barcodeTranslate_hex, cellid)
    
    # Read the beads statistics file and calculate the total number of raw reads and gene reads.
    raw_beads_stat = pd.read_table(beads_stat, sep='\t')
    Total_reads = raw_beads_stat['Raw'].sum()
    Total_gn_reads = raw_beads_stat['GnReads'].sum()
    
    # Read the barcode-to-cell ID translation file and merge it with the beads statistics file.
    barcodeTranslate = pd.read_table(barcodeTranslate, sep='\t', header=None)
    barcodeTranslate.columns = ['BARCODE', 'CELL']
    beads_stat = pd.merge(raw_beads_stat, barcodeTranslate, how='inner', on='BARCODE')
    
    # Drop the BARCODE and GN columns from the merged dataframe, 
    # group the data by cell ID, and aggregate the values by sum.
    beads_stat = beads_stat.drop(['BARCODE', 'GN'], axis=1)
    beads_stat = beads_stat.groupby("CELL").agg('sum')
    
    # Calculate the total number of reads and gene reads per cell, 
    # as well as the number of cells and the mean reads per cell.
    cell_reads = beads_stat['Raw'].sum()
    cell_gn_reads = beads_stat['GnReads'].sum()
    cell_number = beads_stat.shape[0]
    cell_mean_reads = str(round(cell_reads / cell_number))
    
    # Calculate the fraction of gene reads in cells as a percentage.
    Fraction_Reads_ratio = str(round(int(cell_gn_reads) * 100 / int(Total_gn_reads), 2)) + '%'
    
    # Write the summary report to the cellCount file.
    count_report = open(cellCount, 'w')
    count_report.write('Fraction Reads in Cells,%s' % Fraction_Reads_ratio + '\n')
    count_report.write('Estimated Number of Cells,%s' % str(cell_number) + '\n')
    count_report.write('Total Reads Number of Cells,%s' % str(cell_reads) + '\n')
    count_report.write('Mean reads per cell,%s' % str(cell_mean_reads) + '\n')
    count_report.close()

    # Extract the cell number frequencies from the barcode-to-cell ID translation file, 
    # and calculate the number of cells with each frequency.
    barcodeTranslate['frequence'] = barcodeTranslate['CELL'].str.split('_N', expand=True)[1]
    figtable = barcodeTranslate.frequence.value_counts()
    figtable = figtable.reset_index(level=None, drop=False, name=None, inplace=False)
    figtable['index'] = figtable['index'].astype(int)
    figtable['Count'] = figtable.apply(lambda x: round(x['frequence'] / x['index']), axis=1)
    figtable.columns = ['Num', 'frequence', 'Count']
    figtable['Num'] = figtable['Num'].astype(str)
    cellnum = figtable['Count'].sum()
    figtable['num_count'] = figtable["Num"].map(str) +'  '+figtable["Count"].map(str)
    figtable = figtable.sort_values("Num")
    merge_graph(figtable,cellnum,outdir)

def parse_args():
    parser = argparse.ArgumentParser(description='summary barcode and cell merge')
    parser.add_argument('--combined_file', type=str, help='combined_list file for analysis')
    parser.add_argument('--select_barcode', type=str, help='select barcode')
    parser.add_argument('--beads_stat', type=str, help='beads stat')
    parser.add_argument('--outdir', type=str, help='set the outdir for analysis')
    parser.add_argument('--name',type=str, help='sample name')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    combined_file = args.combined_file
    select_barcode = args.select_barcode
    beads_stat = args.beads_stat
    outdir = args.outdir
    name= args.name
    barcodeTranslate = '%s/%s_barcodeTranslate.txt'%(outdir,name)
    barcodeTranslate_hex = '%s/%s_barcodeTranslate_hex.txt'%(outdir,name)
    cellCount = '%s/cellCount_report.csv'%outdir
    cellid = '%s/cell.id'%outdir
    summary_count(combined_file,select_barcode,beads_stat,barcodeTranslate,barcodeTranslate_hex,cellid,cellCount,outdir)

if __name__ == '__main__':
    main()