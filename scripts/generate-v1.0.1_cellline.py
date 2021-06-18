import pandas as pd
import os
import argparse
import re

def is_number(s):
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
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    parser.add_argument('--htmlTemplate', type=str, help=
	'''input the html template''',)
    parser.add_argument('--ID', type=str, help=
        '''input the ID''',)
    args = parser.parse_args()
    return [args.outPath, args.htmlTemplate, args.ID,]
 
def get_args_from_file():
    path=get_args()[0]
    csv = [path+'/07.report/1.cell_report.csv',\
    path+'/07.report/3.cDNA_sequencing_report.csv',\
    path+'/07.report/3.Index_sequencing_report_T1.csv',\
    path+'/07.report/4.alignment_report.csv',\
    path+'/07.report/5.anno_report.csv',\
    ]
    
    stat = dict()
    for i in range(len(csv)):
        if i==0:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['samplename'] = df[1][0]
            stat['species'] = df[1][1]
            stat['estm_Num_cell'] = df[1][2]
            stat['mean_r_per_c'] = df[1][3]
            stat['mean_UMI_per_c'] = df[1][4]
            stat['median_UMI_per_c'] = df[1][5]
#            stat['total_gene'] = df[1][6]
            
            stat['mean_genes_per_c'] = df[1][6]
            stat['median_genes_per_c'] = df[1][7]
            stat['Human_mean_genes'] = df[1][8]
            stat['Mouse_mean_genes'] = df[1][9]          
            stat['Human_mean_UMIs'] = df[1][10]            
            stat['Mouse_mean_UMIs'] = df[1][11]
        if i==1:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['cDNA_num_frag'] = df[1][0]
            stat['cDNA_frag_pass_QC'] = df[1][1]
            stat['cDNA_frag_exact_bar'] = df[1][2]
            stat['cDNA_frag_fail_bar'] = df[1][3]
            stat['cDNA_frag_low_qual'] = df[1][4]
            stat['cDNA_frag_unknow_bar'] = df[1][5]
            stat['cDNA_Q30_c_bar'] = df[1][6]
            stat['cDNA_Q30_s_bar'] = df[1][7]
            stat['cDNA_Q30_UMI'] = df[1][8]
            stat['cDNA_Q30_r'] = df[1][9]
        if i==2:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['index_num_frag'] = df[1][0]
            stat['index_frag_pass_QC'] = df[1][1]
            stat['index_frag_exact_bar'] = df[1][2]
            stat['index_frag_fail_bar'] = df[1][3]
            stat['index_frag_low_qual'] = df[1][4]
            stat['index_frag_unknow_bar'] = df[1][5]
            stat['index_Q30_c_bar'] = df[1][6]
            stat['index_Q30_s_bar'] = df[1][7]
            stat['index_Q30_UMI'] = df[1][8]
            stat['index_Q30_r'] = df[1][9]
        if i==3:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")        
            stat['raw_r'] = df[1][0]
            stat['map_r'] = df[1][1]
            stat['plus_strd'] = df[1][2]
            stat['minus_strd'] = df[1][3]
            stat['mito_ratio'] = df[1][4]
            stat['map_qual_corrt_r'] = df[1][5]            
        if i==4:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['r_m_geno'] = df[1][0]
            stat['r_m_ex'] = df[1][1]
            stat['r_m_intro'] = df[1][2]
            stat['r_m_ex_intro'] = df[1][3]
            stat['r_m_anti'] = df[1][4]
            stat['r_m_inter'] = df[1][5]
            stat['r_m_gene_fail'] = df[1][6]
    plot_file = [
    path+'/07.report/div/barcode_rank.div',\
    path+'/07.report/div/cluster_chsize.div',\
    path+'/07.report/base64/6.base64',\
    path+'/07.report/base64/7.base64',\
    #path+'/07.report/base64/8.base64',\
    
]
    plot_base64 = []
    plot_base64.append(open(path+'/07.report/div/barcode_rank.div',"r").read())
#    plot_base64.append(open(path+'/07.report/div/cluster_chsize.div',"r").read())
    plot_base64.append(open(path+'/07.report/base64/6.base64',"r").read())
    plot_base64.append(open(path+'/07.report/base64/7.base64',"r").read())
    #plot_base64.append(open(path+'/07.report/base64/8.base64',"r").read())
#    plot_base64.append(open(path+'/07.report/div/nUMI_chsize.div',"r").read())
    plot_order = ['plot1','plot3','plot4']
    plot_dict = dict(zip(plot_order, plot_base64))
    
    import locale
    locale.setlocale(locale.LC_ALL, 'en_US')
    for k,v in stat.items():
        if is_number(v):
            stat[k] =locale.format_string("%d", int(v), grouping=True)
        else:
            continue
    return stat, plot_dict
    
def write_param_to_template():
    stat, plot_dict= get_args_from_file()
    template = open(get_args()[1]).read()
    ID = get_args()[2]
    from string import Template

    html=Template(template)
    report=html.safe_substitute(sample_info=ID, samplename=stat['samplename'],\
                    species=stat['species'], median_UMI_per_c=stat['median_UMI_per_c'],\
                    estm_Num_cell=stat['estm_Num_cell'],sample_id=stat['estm_Num_cell'],\
                  #  total_gene=stat['total_gene'],\
                  #  cluster_cell=stat['cluster_cell'],\
                    cDNA_num_frag=stat['cDNA_num_frag'],\
                    cDNA_frag_pass_QC=stat['cDNA_frag_pass_QC'],cDNA_frag_exact_bar=stat['cDNA_frag_exact_bar'],\
                    cDNA_frag_fail_bar=stat['cDNA_frag_fail_bar'],\
                    cDNA_frag_low_qual=stat['cDNA_frag_low_qual'],\
                    cDNA_frag_unknow_bar = stat['cDNA_frag_unknow_bar'],\
                    cDNA_Q30_c_bar=stat['cDNA_Q30_c_bar'],cDNA_Q30_s_bar=stat['cDNA_Q30_s_bar'],\
                    cDNA_Q30_UMI=stat['cDNA_Q30_UMI'],\
                    cDNA_Q30_r=stat['cDNA_Q30_r'],\
                    index_num_frag=stat['index_num_frag'],\
                    index_frag_pass_QC=stat['index_frag_pass_QC'],index_frag_exact_bar=stat['index_frag_exact_bar'],\
                    index_frag_fail_bar=stat['index_frag_fail_bar'],\
                    index_frag_low_qual=stat['index_frag_low_qual'],\
                    index_frag_unknow_bar = stat['index_frag_unknow_bar'],\
                    index_Q30_c_bar=stat['index_Q30_c_bar'],index_Q30_s_bar=stat['index_Q30_s_bar'],\
                    index_Q30_UMI=stat['index_Q30_UMI'],\
                    index_Q30_r=stat['index_Q30_r'],\
                    raw_r=stat['raw_r'],\
                    map_r=stat['map_r'],plus_strd=stat['plus_strd'],\
                    minus_strd=stat['minus_strd'],\
                    mito_ratio = stat['mito_ratio'],\
                    map_qual_corrt_r=stat['map_qual_corrt_r'],plot1=plot_dict['plot1'],\
                    plot3=plot_dict['plot3'],\
                    plot4=plot_dict['plot4'],r_m_geno=stat['r_m_geno'], 
            r_m_ex=stat['r_m_ex'], 
            r_m_intro=stat['r_m_intro'],
            r_m_ex_intro=stat['r_m_ex_intro'],
            r_m_anti=stat['r_m_anti'],
            r_m_inter=stat['r_m_inter'],
            r_m_gene_fail=stat['r_m_gene_fail'],ID=ID,
            mean_r_per_c=stat['mean_r_per_c'],
            mean_UMI_per_c=stat['mean_UMI_per_c'],
            mean_genes_per_c=stat['mean_genes_per_c'],
            median_genes_per_c=stat['median_genes_per_c'],
            #plot6=plot_dict['plot6']
	            
            #mean_genes_per_c=stat['mean_genes_per_c'],
            #median_genes_per_c=stat['median_genes_per_c'],
            Human_mean_genes=stat['Human_mean_genes'],
            Mouse_mean_genes=stat['Mouse_mean_genes'],         
            Human_mean_UMIs=stat['Human_mean_UMIs'],           
            Mouse_mean_UMIs=stat['Mouse_mean_UMIs'],

                    )
    return report
    
if __name__ == '__main__':
    outpath=get_args()[0]
    ID = get_args()[2]
    get_args_from_file()
    report=write_param_to_template()
    fw = open(outpath+'/07.report/html/'+ID+'_CDCPv2_scRNA_report.html','w')
    fw.write(report)
