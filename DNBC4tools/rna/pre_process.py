import plotly.express as px
import pandas as pd
import plotly as py
import plotly.graph_objs as go
from plotly.io import *
import os,collections,argparse,math
import numpy as np
from scipy.interpolate import make_interp_spline

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    parser.add_argument('--sample', type=str, help=
	'''input the samplename''',)
    args = parser.parse_args()
    return args.outPath,args.sample
    
outpath,samplename = get_args()

os.system('mkdir -p %s' % (outpath+"/04.report/div"))
os.system('mkdir -p %s' % (outpath+"/04.report/base64"))
os.system('mkdir -p %s' % (outpath+"/04.report/table"))
#os.system('mkdir -p %s' % (outpath+"/04.report/html"))

config = {
    'modeBarButtonsToRemove': ["autoScale2d","hoverClosestCartesian", "hoverCompareCartesian", "lasso2d",
    "zoomIn2d", "zoomOut2d", "sendDataToCloud",
    "toggleSpikelines" ,"logo"],
    'displaylogo': False,}

print('''############################
# 1.barcode plot
############################''')
config=config

def segment_log_plot(y_data, x_start, x_end):
    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))
    segment_len = 0.0
    segment_idx = [x_start]
    for i in range(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        segment_len += np.linalg.norm([dx, dy])
        if segment_len >= 0.02 and i > (segment_idx[-1] + 20):
            segment_idx.append(i+1)
            segment_len = 0.0
    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)
    return segment_idx

def plot_cmap(density):
    plot_colors =  ["#DDDDDD","#D6D9DC","#CFD6DB","#C8D3DA","#C1D0D9","#BACDD9","#B3C9D8","#ACC6D7","#A5C3D6","#9EC0D6",
                    "#97BDD5","#90BAD4","#89B6D3","#82B3D3","#7BB0D2","#74ADD1","#6DAAD0","#66A6CF","#5FA3CF","#58A0CE",
                    "#539DCC","#4F99CA","#4C95C8","#4992C6","#458EC3","#428AC1","#3F87BF","#3B83BD","#3880BA","#347CB8",
                    "#3178B6","#2E75B4","#2A71B1","#276DAF","#236AAD","#2066AB","#1D62A8","#195FA6","#165BA4","#1358A2"]
    levels = len(plot_colors)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return plot_colors[ind]

def plot_jaccard_knee_frag(count_data_path):
    count_data = pd.read_csv(count_data_path, index_col=0, sep=',')
    count_data['New']=count_data.index
    cell_bc = np.array(count_data[count_data['Beads'] == 'true'].index)
    sorted_bc = np.array(count_data.index)
    sorted_counts = np.array(count_data['UMI'])
    total_bc = len(sorted_bc)
    ix1 = count_data.drop_duplicates('Beads',keep='first').index[1]-1
    ix2 = count_data.drop_duplicates('Beads',keep='last').index[0]
    plot_segments = []
    barcodeSegment = collections.namedtuple('barcodeSegment', ['start', 'end', 'density', 'legend'])
    plot_segments.append(barcodeSegment(start=0, end=ix1, density=1.0, legend=True))
    plot_segments.append(barcodeSegment(start=ix2+1, end=total_bc, density=0.0, legend=True))
    mixed_segments = segment_log_plot(sorted_counts, ix1, ix2)
    for i in range(len(mixed_segments) - 1):
        num_cells = sum([1 for i in range(mixed_segments[i], mixed_segments[i + 1]) if sorted_bc[i] in cell_bc])
        density = float(num_cells)/float(mixed_segments[i + 1]-mixed_segments[i])
        plot_segments.append(barcodeSegment(start=mixed_segments[i], end = mixed_segments[i + 1], density=density, legend=False))
    plot_data = []
    for plot_segment in plot_segments:
        start = max(0, plot_segment.start - 1)
        end = plot_segment.end
        selct_count = count_data[start:end]
        dp_first = set(selct_count[selct_count[["UMI"]].duplicated(keep="first")].index)
        dp_last = set(selct_count[selct_count[["UMI"]].duplicated(keep="last")].index)
        dp_inter = dp_first & dp_last
        selct_count=selct_count.drop(list(dp_inter),axis=0)
        x = list(selct_count['New'])
        y = list(selct_count['UMI'])
        name = 'TRUE' if plot_segment.density > 0 else 'NOISE'
        if plot_segment.density > 0:
            n_barcodes = plot_segment.end - plot_segment.start
            n_cells = int(round(plot_segment.density * n_barcodes))
            hover = "{:.0f}% Beads<br>({}/{})".format(100 * plot_segment.density, n_cells, n_barcodes)
        else:
            hover = "NOISE"
        data_dict = {"x": x,"y": y,"name": name,"hoverinfo": "text","text": hover,"type": "scattergl","mode": "lines","line": {"width": 3,"color": plot_cmap(plot_segment.density),},"showlegend": plot_segment.legend,}
        plot_data.append(data_dict)
    plotly_data = [go.Scatter(x=dat['x'], y=dat['y'], name=dat['name'], mode=dat['mode'], showlegend=dat['showlegend'],marker={'color': dat['line']['color']}, line=dat['line'], text=dat['text']) for dat in plot_data]
    layout = go.Layout(xaxis=dict(type="log", gridcolor="lightgrey",title="Barcode in Rank-descending Order",color="black",showline=True,zeroline=True,linewidth=1,fixedrange= True,linecolor="black"),
                    yaxis = dict(type="log",title="UMI counts",gridcolor="lightgrey",linewidth=1,fixedrange= True,color="black",linecolor="black"),
                    height=360,width=450,plot_bgcolor='rgba(0,0,0,0)',hovermode='closest',paper_bgcolor='white',legend=dict(x=1,y=1,traceorder="normal",font=dict(family="Arial",size=12,color="black"),bordercolor="Black",borderwidth=0),
                    margin=dict(l=0,r=0,b=0,t=0,pad=1),font=dict(size=10))
    fig = go.Figure(data=plotly_data, layout=layout)
    py.offline.plot(fig, filename = outpath+"/04.report/div/barcode_rank.html",auto_open=False,config=config)
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(outpath+"/04.report/div/barcode_rank.div",'w')
    fw.write(fig2)
    fw.close()

plot_jaccard_knee_frag(outpath+"/02.count/cutoff.csv")

print('''############################
# 2.cluster plot
############################''')
outdir=outpath
cluster_file = outdir+'/03.analysis/Clustering/cluster.csv'
cluster = pd.read_csv(cluster_file)
df = pd.read_csv(cluster_file)
#df = df.sort_values(by='Cluster')
df[['Cluster']] = df[['Cluster']].astype('str')    
fig = px.scatter(df, x=df.UMAP_1, y=df.UMAP_2, color= df['Cluster'],color_discrete_sequence=px.colors.qualitative.G10)

config=config
fig.update_layout(
    autosize=False,
    width=565,
    height=500,
    legend_title=dict(font=dict(size=13),text='Cluster',),
    legend=dict(font=dict(size=10,family='Arial'),itemsizing='constant'),
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',),
    yaxis=dict(gridcolor='lightgray',)
    )
fig.update_traces(marker={'size': 3})
fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray') 
py.offline.plot(fig, filename = outdir+"/04.report/div/cluster.html",auto_open=False,config=config)
fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
fw = open(outdir+"/04.report/div/cluster.div",'w')
fw.write(fig2)
fw.close()
barplot_chsize = open(outdir+"/04.report/div/cluster.div","r").read()
barplot_chsize=barplot_chsize.replace('width:100%','width:565px')
fw1 = open(outdir+"/04.report/div/cluster_chsize.div","w")
fw1.write(barplot_chsize)
fw1.close()

####################plot nUMI
fig = px.scatter(df, x=df.UMAP_1, y=df.UMAP_2, color= df['nUMI'], )
config=config
fig.update_layout(
    autosize=False,
    width=520,
    height=500,
    legend_title=dict(font=dict(size=20),text='nUMI',),
    legend=dict(font=dict(size=20,family='Arial'),),
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',),
    yaxis=dict(gridcolor='lightgray',)
    )
fig.update_traces(marker={'size': 3})
fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
py.offline.plot(fig, filename = outdir+"/04.report/div/cluster_nUMI.html",auto_open=False,config=config)
fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
fw = open(outdir+"/04.report/div/nUMI.div",'w')
fw.write(fig2)
fw.close()
barplot_chsize = open(outdir+"/04.report/div/nUMI.div","r").read()
barplot_chsize=barplot_chsize.replace('width:100%','width:480px')
fw1 = open(outdir+"/04.report/div/nUMI_chsize.div","w")
fw1.write(barplot_chsize)
fw1.close()


print('''###########################
# 2.2 cluster annotation plot
###########################''')
####################plot nUMI
if 'Predicted cell type' in df.columns:
    fig = px.scatter(df, x=df.UMAP_1, y=df.UMAP_2, color= df['Predicted cell type'],color_discrete_sequence=px.colors.qualitative.G10)
    config=config
    fig.update_layout(
    autosize=False,
    width=900,
    height=500,
    legend_title=dict(font=dict(size=16),text='Predicted cell type: cell number',),
    legend=dict(x=1.2,y=0.5,font=dict(size=10,family='Arial'),itemsizing='constant'),
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',),
    yaxis=dict(gridcolor='lightgray',)
    )
    fig.update_traces(marker={'size': 3})
    fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
    fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
    py.offline.plot(fig, filename = outdir+"/04.report/div/anno.html",auto_open=False,config=config)
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(outdir+"/04.report/div/anno.div",'w')
    fw.write(fig2)
    fw.close()
    barplot_chsize = open(outdir+"/04.report/div/anno.div","r").read()
    barplot_chsize=barplot_chsize.replace('width:100%','width:565px')
    fw1 = open(outdir+"/04.report/div/anno_chsize.div","w")
    fw1.write(barplot_chsize)
    fw1.close()

print('''###########################
# 2.3 saturation plot
###########################''')
####################plot nUMI
outdir=outpath
saturation_file = outdir+'/02.count/'+'saturation.xls'
df = pd.read_csv(saturation_file, sep="\t",header=0)

x=df['Mean Reads per Cell']
y=df['Sequencing Saturation']
if len(df) > 2:
    xnew = np.linspace(x.min(),x.max(),300)
    #import statsmodels.api as sm
    #lowess = sm.nonparametric.lowess
    #ynew = lowess(y, x, frac=0.27)
    ynew = make_interp_spline(x,y)(xnew)
    #fig = px.line(df, x=ynew[:,0], y=ynew[:,1])
else:
    xnew = x
    ynew = y
fig = px.line(df, x=xnew, y=ynew )

config=config
fig.update_layout(
    autosize=False,
    width=565,
    height=500,
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',title="Mean Reads per Cell"),
    yaxis=dict(gridcolor='lightgray',title="Sequencing Saturation"),
    yaxis_range=[0,100]
    )
fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_traces(line=dict(color="#337ab7", width=3))
py.offline.plot(fig, filename = outdir+"/04.report/div/saturation.html",auto_open=False,config=config)
fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
fw = open(outdir+"/04.report/div/saturation.div",'w')
fw.write(fig2)
fw.close()

x=df['Mean Reads per Cell']
y=df['Median Genes per Cell']
if len(df) > 2:
    xnew = np.linspace(x.min(),x.max(),400)
    ynew = make_interp_spline(x,y)(xnew)
    #ynew = lowess(y, x, frac=0.27)
else:
    xnew = x
    ynew = y
#fig = px.line(df, x=ynew[:,0], y=ynew[:,1])
fig = px.line(df, x=xnew, y=ynew )
config=config
fig.update_layout(
    autosize=False,
    width=520,
    height=500,
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',title="Mean Reads per Cell"),
    yaxis=dict(gridcolor='lightgray',title="Median Genes per Cell")
    )
fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_traces(line=dict(color="#337ab7", width=3))
py.offline.plot(fig, filename = outdir+"/04.report/div/saturation2.html",auto_open=False,config=config)
fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
fw = open(outdir+"/04.report/div/saturation2.div",'w')
fw.write(fig2)
fw.close()

print('''###########################
# 3.png to base64
###########################''')
import base64
def png_to_base64(file=str(),filename=str()):    
    #inpath = get_args() 
    outpath = outdir + "/04.report"
    file_path = outdir+"/"+file
    base64_path = outpath+'/base64'+'/'+filename+'.base64'
    if os.path.isfile(file_path):
        with open(file_path, "rb") as f:
            base64_data = base64.b64encode(f.read())
            s = base64_data.decode()
            base64_path_f = open(base64_path, 'w')
            base64_path_f.write('<img src=data:image/'+'png'+';base64,'+s+">")
            base64_path_f.close()           
pictures = {'02.count/cellNumber_merge.png':'6','03.analysis/QC/raw_QCplot.png':'7'}
for k,v in pictures.items():
    png_to_base64(k,v) 

print('''###########################
# 4.csv to data-table
###########################''')
if os.path.exists(outdir+'/03.analysis/Clustering/marker.csv'): 
    df1= pd.read_csv(open(outdir+'/03.analysis/Clustering/marker.csv'),encoding="utf-8",dtype=str,)
    fw = open(outdir+'/04.report/table/marker-table.txt','w')
    for index, row in df1.iterrows():
        fw.write('<tr><td>'+row['gene']+'</td>'\
                +'<td>'+row['cluster']+'</td>'\
                +'<td>'+row['p_val_adj']+'</td>'\
                +'<td>'+row['p_val']+'</td>'\
                +'<td>'+row['avg_log2FC']+'</td>'\
                +'<td>'+row['pct.1']+'</td>'\
                +'<td>'+row['pct.2']+'</td>'\
            )
