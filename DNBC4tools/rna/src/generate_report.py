import pandas as pd
import os,argparse
from dnbc4tools.__init__ import __version__

def is_number(s):
    """
    Check if the given string can be converted to a number.
    :param s: A string to check.
    :return: True if the string can be converted to a number, False otherwise.
    """
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False
 
def get_args_from_file(path,intron,species,samplename):
    csv = [
        path+'/02.count/cellCount_report.csv',
        path+'/02.count/saturation.xls',
        path+'/01.data/cDNA_sequencing_report.csv',
        path+'/01.data/Index_sequencing_report.csv',
        path+'/01.data/alignment_report.csv',
        path+'/01.data/anno_report.csv',
        ]
    stat = dict()
    stat['intron_boolean'] = intron
    stat['samplename'] = samplename
    stat['species'] = species
    stat['version'] = __version__
    for i in range(len(csv)):
        if i==0:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['Fraction'] = df[1][0]
            stat['estm_Num_cell'] = df[1][1]
            stat['mean_r_per_c'] = df[1][3]
            stat['mean_UMI_per_c'] = df[1][4]
            stat['median_UMI_per_c'] = df[1][5]
            stat['total_gene'] = df[1][6]
            stat['mean_genes_per_c'] = df[1][7]
            stat['median_genes_per_c'] = df[1][8]
        if i==1:
            df = pd.read_table(open(csv[i]),encoding="utf_8",dtype=str,header=0,sep="\t",index_col=0)
            stat['saturation'] = str(df.iloc[[11],[4]].values[0][0])+'%'
        if i==2:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['cDNA_num_frag'] = df[1][0]
            stat['cDNA_frag_pass_QC'] = df[1][1]
            stat['cDNA_frag_low_qual'] = df[1][2]
            stat['cDNA_frag_fail_bar'] = df[1][3]
            stat['cDNA_frag_exact_bar'] = df[1][5]
            stat['cDNA_adapter'] = df[1][6]
            stat['cDNA_Q30_c_bar'] = df[1][7]
            stat['cDNA_Q30_s_bar'] = df[1][8]
            stat['cDNA_Q30_UMI'] = df[1][9]
            stat['cDNA_Q30_r'] = df[1][10]
        if i==3:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['index_num_frag'] = df[1][0]
            stat['index_frag_pass_QC'] = df[1][1]
            stat['index_frag_low_qual'] = df[1][2]
            stat['index_frag_fail_bar'] = df[1][3]
            stat['index_frag_exact_bar'] = df[1][5]
            stat['index_Q30_c_bar'] = df[1][7]
            stat['index_Q30_s_bar'] = df[1][8]
            stat['index_Q30_UMI'] = df[1][9]
            stat['index_Q30_r'] = df[1][10]
        if i==4:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")        
            stat['raw_r'] = df[1][0]
            stat['map_r'] = df[1][1]
            stat['plus_strd'] = df[1][2]
            stat['minus_strd'] = df[1][3]
            stat['mito_ratio'] = df[1][4]
            stat['map_qual_corrt_r'] = df[1][5]         
        if i==5:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['r_m_ex'] = df[1][1]
            stat['r_m_intro'] = df[1][2]
            stat['r_m_ex_intro'] = df[1][3]
            stat['r_m_anti'] = df[1][4]
            stat['r_m_inter'] = df[1][5]
            stat['r_m_gene_fail'] = df[1][6]
            stat['r_m_geno'] = str(round(int(stat['raw_r'])*100/int(stat['cDNA_frag_pass_QC']),2))+'%'

    stat['cDNA_frag_pass_QC_ratio'] = str(round(int(stat['cDNA_frag_pass_QC'])*100/int(stat['cDNA_num_frag']),2))+'%'
    stat['cDNA_frag_exact_bar'] = str(round(int(stat['cDNA_frag_exact_bar'])*100/int(stat['cDNA_num_frag']),2))+'%'
    stat['cDNA_frag_fail_bar'] = str(round(int(stat['cDNA_frag_fail_bar'])*100/int(stat['cDNA_num_frag']),2))+'%'
    stat['cDNA_frag_low_qual'] = str(round(int(stat['cDNA_frag_low_qual'])*100/int(stat['cDNA_num_frag']),2))+'%'
    stat['index_frag_pass_QC_ratio'] = str(round(int(stat['index_frag_pass_QC'])*100/int(stat['index_num_frag']),2))+'%'
    stat['index_frag_exact_bar'] = str(round(int(stat['index_frag_exact_bar'])*100/int(stat['index_num_frag']),2))+'%'
    stat['index_frag_fail_bar'] = str(round(int(stat['index_frag_fail_bar'])*100/int(stat['index_num_frag']),2))+'%'
    stat['index_frag_low_qual'] = str(round(int(stat['index_frag_low_qual'])*100/int(stat['index_num_frag']),2))+'%'
    stat['plus_strd'] = str(round(int(stat['plus_strd'])*100/int(stat['raw_r']),2))+'%'
    stat['minus_strd'] = str(round(int(stat['minus_strd'])*100/int(stat['raw_r']),2))+'%'
    stat['map_qual_corrt_r'] = str(round(int(stat['map_qual_corrt_r'])*100/int(stat['raw_r']),2))+'%'

    plot_file = [
        path+'/04.report/div/barcode_rank.div',
        path+'/04.report/div/cluster_chsize.div',
        path+'/04.report/base64/6.base64',
        path+'/04.report/base64/7.base64',
        path+'/04.report/div/nUMI_chsize.div',
        path+'/04.report/div/saturation.div',
        path+'/04.report/div/saturation2.div',
        path+'/04.report/div/anno_chsize.div',
        ]

    plot_base64 = []
    plot_base64.append(open(path+'/04.report/div/barcode_rank.div',"r").read())
    plot_base64.append(open(path+'/04.report/div/cluster_chsize.div',"r").read())
    plot_base64.append(open(path+'/04.report/base64/6.base64',"r").read())
    plot_base64.append(open(path+'/04.report/base64/7.base64',"r").read())
    plot_base64.append(open(path+'/04.report/div/nUMI_chsize.div',"r").read())
    plot_base64.append(open(path+'/04.report/div/saturation.div',"r").read())
    plot_base64.append(open(path+'/04.report/div/saturation2.div',"r").read())

    if os.path.exists(path+'/04.report/div/anno_chsize.div'):
        plot_base64.append(open(path+'/04.report/div/anno_chsize.div',"r").read())
    
    '''
    for f in plot_file:
        if re.search('6.base64',f) or re.search('7.base64',f) :
            if os.path.exists(f):
                base64 = open(f).read()
                #img = ("<img src=%s height=500px width=100\%>" %base64)
                img = "<img src=\"data:image/png+xml;base64,"+base64+"\">"
                plot_base64.append(img)
            else:
                plot_base64.append(
                <p style="font_family=DIN Next LT Pro;font_size=18px;font_weight=400">
                The cluster plot has not been generated because the data quality is too low.
                <p>
                )
        else:
            plot_base64.append(open(f).read())
    '''
    
    plot_order = ['plot1','plot2','plot3','plot4','plot5','plot6','plot7', 'plot8']
    plot_dict = dict(zip(plot_order, plot_base64))
     
    data_tables_file = path+'/04.report/table/marker-table.txt'
    table = str()
    if os.path.exists(data_tables_file):
        table='''<table id=\"table_id_example\" class=\"table table-bordered table-striped\" style=\"Dosis;\">            <thead style=\"font-size:11px\"><tr>
                <th>gene</th>
                <th>cluster</th>
                <th>p_val_adj</th>
                <th>p_val</th>
                <th>avg_log2FC</th>
                <th>pct.1</th>
                <th>pct.2</th>
                

                

            </tr>
        </thead>
            <tbody style=\"font-size:11px;\">
            '''+open(data_tables_file).read()+"</tbody></table>"
            #data_tables.append(table)               
        
    else:
        table.append('''
        <p style="font_family=DIN Next LT Pro;font_size=18px;font_weight=400">
        The table has not been generated because the data quality is too low.
        <p>
        ''')

    # import locale
    # locale.setlocale(locale.LC_ALL, 'en_US')
    for k,v in stat.items():
        if is_number(v):
            # stat[k] =locale.format_string("%d", int(v), grouping=True)
            stat[k] = format(int(v),',')
        else:
            continue
    return stat, plot_dict, table
    
def write_param_to_template(htmlTemplate,samplename,path,intron,species):
    stat, plot_dict, table = get_args_from_file(path,intron,species,samplename)
    template = open(htmlTemplate).read()
    from string import Template
    html=Template(template)
    if os.path.exists(path+'/04.report/div/anno_chsize.div'):
        plot8_content = plot_dict['plot8']
    else:
        plot8_content = "There is no such species reference for annnotation."

    report=html.safe_substitute(
            sample_info=samplename, 
            samplename=stat['samplename'],
            species=stat['species'], 
            median_UMI_per_c=stat['median_UMI_per_c'],
            estm_Num_cell=stat['estm_Num_cell'],
            sample_id=stat['estm_Num_cell'],
            total_gene=stat['total_gene'], 
            saturation=stat['saturation'],
            cDNA_num_frag=stat['cDNA_num_frag'],
            intron_boolean=stat['intron_boolean'],
            cDNA_frag_pass_QC=stat['cDNA_frag_pass_QC_ratio'],
            cDNA_frag_exact_bar=stat['cDNA_frag_exact_bar'],
            cDNA_frag_fail_bar=stat['cDNA_frag_fail_bar'],
            cDNA_frag_low_qual=stat['cDNA_frag_low_qual'],
            cDNA_Q30_c_bar=stat['cDNA_Q30_c_bar'],
            cDNA_Q30_s_bar=stat['cDNA_Q30_s_bar'],
            cDNA_Q30_UMI=stat['cDNA_Q30_UMI'],
            cDNA_Q30_r=stat['cDNA_Q30_r'],
            index_num_frag=stat['index_num_frag'],
            Fraction=stat['Fraction'],
            index_frag_pass_QC=stat['index_frag_pass_QC_ratio'],
            index_frag_exact_bar=stat['index_frag_exact_bar'],
            index_frag_fail_bar=stat['index_frag_fail_bar'],
            index_frag_low_qual=stat['index_frag_low_qual'],
            index_Q30_c_bar=stat['index_Q30_c_bar'],
            index_Q30_s_bar=stat['index_Q30_s_bar'],
            index_Q30_UMI=stat['index_Q30_UMI'],
            index_Q30_r=stat['index_Q30_r'],
            raw_r=stat['cDNA_frag_pass_QC'],
            map_r=stat['raw_r'],
            plus_strd=stat['plus_strd'],
            minus_strd=stat['minus_strd'],
            mito_ratio = stat['mito_ratio'],
            map_qual_corrt_r=stat['map_qual_corrt_r'],
            plot1=plot_dict['plot1'],
            plot2=plot_dict['plot2'],
            plot3=plot_dict['plot3'],
            plot4=plot_dict['plot4'],
            plot5=plot_dict['plot5'],
            table = table, 
            r_m_geno=stat['r_m_geno'],
            plot6=plot_dict['plot6'],
            plot7=plot_dict['plot7'],
            plot8= plot8_content,
            r_m_ex=stat['r_m_ex'],
            r_m_intro=stat['r_m_intro'],
            pipeversion = stat['version'],
            r_m_ex_intro=stat['r_m_ex_intro'],
            r_m_anti=stat['r_m_anti'],
            r_m_inter=stat['r_m_inter'],
            r_m_gene_fail=stat['r_m_gene_fail'],
            ID=samplename,
            mean_r_per_c=stat['mean_r_per_c'],
            mean_UMI_per_c=stat['mean_UMI_per_c'],
            mean_genes_per_c=stat['mean_genes_per_c'],
            median_genes_per_c=stat['median_genes_per_c']
            )
    
    metrics_df = pd.DataFrame([stat])
    cols = ["samplename","species","estm_Num_cell","mean_r_per_c","mean_UMI_per_c","median_UMI_per_c","total_gene",\
        "mean_genes_per_c","median_genes_per_c","saturation","Fraction","cDNA_num_frag","cDNA_frag_pass_QC_ratio","cDNA_adapter","cDNA_Q30_r",\
        "index_num_frag","index_frag_pass_QC_ratio","index_Q30_r","mito_ratio","r_m_geno","r_m_ex","r_m_intro",\
        "r_m_anti","r_m_inter"]
    metrics_summary_df = metrics_df[cols]
    metrics_summary_df.columns =["SampleName","species","Estimated number of cell","Mean reads per cell","Mean UMI count per cell",\
        "Median UMI counts per cell","Total genes detected","Mean genes per cell","Median genes per cell","Sequencing saturation",\
        "Fraction Reads in cell","cDNA Number of reads","cDNA Reads pass QC","cDNA Adapter Reads","cDNA Q30 bases in reads","index Number of reads",\
        "index Reads pass QC","index Q30 bases in reads","Mitochondria ratio","Reads mapped to genome (Map Quality >= 0)",\
        "Reads mapped to exonic regions","Reads mapped to intronic regions","Reads mapped antisense to gene","Reads mapped to intergenic regions"]
    return report,metrics_summary_df

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    parser.add_argument('--htmlTemplate', type=str, help=
	'''input the html template''',)
    parser.add_argument('--name', type=str, help=
    '''input the sample name''',)
    parser.add_argument('--species', type=str, help=
    '''input the sample species, default is undefined.''', default='undefined')
    parser.add_argument('--intron', type=str, help=
    '''True or False of intron reads''')
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    report,metrics_summary_df = write_param_to_template(args.htmlTemplate,args.name,args.outPath,args.intron,args.species)
    fw = open(args.outPath+'/04.report/'+args.name+'_scRNA_report.html','w')
    fw.write(report)
    file_df = args.outPath +'/04.report/metrics_summary.xls'
    metrics_summary_df.to_csv(file_df,sep='\t',index=None)
    
if __name__ == '__main__':
    main()