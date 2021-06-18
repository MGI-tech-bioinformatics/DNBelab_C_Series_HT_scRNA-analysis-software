import plotly.express as px
import pandas as pd
import plotly as py
from plotly.graph_objs import Scatter, Layout, Data, Scattergl
import pandas as pd
import plotly.graph_objs as go
from plotly.io import *
import os
import argparse
import numpy as np
import sys



print('''############################
# 1.barcode plot
############################''')
import plotly as py
from plotly.graph_objs import Scatter, Layout, Data, Scattergl
import pandas as pd
import plotly.graph_objs as go
from plotly.io import *
import os
import argparse
import numpy as np
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    args = parser.parse_args()
    return args.outPath
    
outpath = get_args()

os.system('mkdir -p %s' % (outpath+"/07.report/div"))
os.system('mkdir -p %s' % (outpath+"/07.report/base64"))
os.system('mkdir -p %s' % (outpath+"/07.report/table"))
os.system('mkdir -p %s' % (outpath+"/07.report/html"))

config = {
    'modeBarButtonsToRemove': ["autoScale2d","hoverClosestCartesian", "hoverCompareCartesian", "lasso2d",
    "zoomIn2d", "zoomOut2d", "sendDataToCloud",
    "toggleSpikelines" ,"logo"],
    'displaylogo': False,}
def plot_jaccard_knee_frag():
    outpath = get_args()
    df = pd.read_csv(open(outpath+"/07.report/2.cut.off.csv"),encoding="utf-8",na_filter=False) 
    #print(df)
    dp_first = set(df[df[["UMI"]].duplicated(keep="first")].index)
    dp_last = set(df[df[["UMI"]].duplicated(keep="last")].index)
    dp_inter = dp_first & dp_last
    df=df.drop(list(dp_inter),axis=0)
    
    df_True = df[(df['Beads'] == "true")]
    df_False = df[(df['Beads'] == "noise")]
    df_NA = df[(df['Beads'] == "NA")]
    
    trace0_x = list(df_True['barcodes'])
    trace0_y = list(df_True['UMI'])  

    trace1_x = list(df_False['barcodes'])
    trace1_y = list(df_False['UMI'])
    
    trace2_x = list(df_NA['barcodes'])
    trace2_y = list(df_NA['UMI'])
   
    
    blue_line = list(zip(trace0_x, trace0_y))
    blue_line = [list(i) for i in blue_line]
    
    red_line = list(zip(trace1_x, trace1_y))
    red_line = [list(i) for i in red_line]  
    
    gray_line = list(zip(trace2_x, trace2_y))
    gray_line = [list(i) for i in gray_line]

    trace0 = Scattergl(
        x = trace0_x,
        y = trace0_y,
        mode="lines",
        name="TRUE",
        line=dict(color="#005BAC",width=3) #zi #67267A
    )
   
    
    trace1 = Scattergl(
        x = trace1_x,
        y = trace1_y,
        mode="lines",
        name="NOISE",
        line=dict(color="gray",width=3)
        
    )
    
    trace2 = Scattergl(
        x = trace2_x,
        y = trace2_y,
        mode="lines",
        name="NA",
        line=dict(color="gray",width=3)
        
    )
    
    

    config={
    'modeBarButtonsToRemove': ["autoScale2d","hoverClosestCartesian", "hoverCompareCartesian", "lasso2d",
    "zoomIn2d", "zoomOut2d", "sendDataToCloud",
    "toggleSpikelines" ,"logo"],
    'displaylogo': False,}
    
    data = [trace0, trace1, trace2]
    
    layout = Layout(
                        xaxis=dict(type="log", 
                        gridcolor="lightgrey",
                        title="Barcode in Rank-descending Order",
                        color="black",
                        showline=True,
                        zeroline=True,
                        linewidth=1,fixedrange= True,
                        linecolor="black"
                        ),
                        
                        yaxis = dict(
                        type="log",
                        title="Reads per Barcode",
                        gridcolor="lightgrey",
                        linewidth=1,fixedrange= True,
                        color="black",
                        linecolor="black"
                        ),
                        height=360,width=450,
                        plot_bgcolor='rgba(0,0,0,0)',
                        
                        hovermode='closest',
                        paper_bgcolor='white',
                        
                        legend=dict(
                        x=1,
                        y=1,
                        traceorder="normal",
                        font=dict(
                        family="Arial",
                        size=12,
                        color="black"
                        ),
                        bordercolor="Black",
                        borderwidth=0
                        ),
                        margin=dict(
                        l=0,
                        r=0,
                        b=0,
                        t=0,
                        pad=1
                        ),
                        font=dict(size=10)
    ) 
    
    
    fig = dict(data=data, layout=layout)
    py.offline.plot(fig, filename = outpath+"/07.report/div/barcode_rank.html",auto_open=False,config=config)
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(outpath+"/07.report/div/barcode_rank.div",'w')
    fw.write(fig2)

path = get_args()
plot_jaccard_knee_frag()


print('''############################
# 2.cluster plot
############################''')
'''
outdir=path
cluster_file = outdir+'/07.report/9.cluster.csv'
cluster = pd.read_csv(cluster_file)
df = pd.read_csv(cluster_file)
 
df = df.sort_values(by='Cluster')
 
df[['Cluster']] = df[['Cluster']].astype('str')    
fig = px.scatter(df, x=df.UMAP_1, y=df.UMAP_2, color= df['Cluster'], )

config=config
fig.update_layout(
    autosize=False,
    width=565,
    height=500,
    legend_title=dict(font=dict(size=16),text='Cluster',),
    legend=dict(font=dict(size=10,family='Arial'),),
    #paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',),
    yaxis=dict(gridcolor='lightgray',)


    )
fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
py.offline.plot(fig, filename = outdir+"/07.report/div/cluster.html",auto_open=False,config=config)
fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
fw = open(outdir+"/07.report/div/cluster.div",'w')
fw.write(fig2)
fw.close()
barplot_chsize = open(outdir+"/07.report/div/cluster.div","r").read()
#print(barplot_chsize)
barplot_chsize=barplot_chsize.replace('width:100%','width:565px')
#print(barplot_chsize)
fw1 = open(outdir+"/07.report/div/cluster_chsize.div","w")
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
    #paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='lightgray',),
    yaxis=dict(gridcolor='lightgray',)


    )
fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
py.offline.plot(fig, filename = outdir+"/07.report/div/cluster_nUMI.html",auto_open=False,config=config)
fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
fw = open(outdir+"/07.report/div/nUMI.div",'w')
fw.write(fig2)
fw.close()
barplot_chsize = open(outdir+"/07.report/div/nUMI.div","r").read()
#print(barplot_chsize)
barplot_chsize=barplot_chsize.replace('width:100%','width:480px')
#print(barplot_chsize)
fw1 = open(outdir+"/07.report/div/nUMI_chsize.div","w")
fw1.write(barplot_chsize)
fw1.close()
'''
print('''###########################
# 3.png to base64
###########################''')
import base64
def png_to_base64(file=str()):    
    outpath = get_args() 
    outpath = outpath + "/07.report"
    filename = file.split('.')[0]
    format = file.split('.')[-1]
    file_path = outpath+"/"+file
    base64_path = outpath+'/base64'+'/'+filename+'.base64'
    if os.path.isfile(file_path):
        with open(file_path, "rb") as f:
            base64_data = base64.b64encode(f.read())
            s = base64_data.decode()
            base64_path_f = open(base64_path, 'w')
            if format =="png":
                base64_path_f.write('<img src=data:image/'+format+';base64,'+s+">")           
                #base64_path_f.write(s) 
                
pictures = [] 
for root,dir,files in os.walk(path):
    for f in files:
        if f.split('.')[-1] == "png":
            pictures.append(f)
 
for p in pictures:
    png_to_base64(p) 

