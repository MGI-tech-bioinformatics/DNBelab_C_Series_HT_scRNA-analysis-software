import plotly.express as px
import pandas as pd
import plotly as py
import plotly.graph_objs as go
from plotly.io import *
import collections,math
import numpy as np
from scipy.interpolate import make_interp_spline

config = {
    'modeBarButtonsToRemove': [
        "autoScale2d",
        "hoverClosestCartesian", 
        "hoverCompareCartesian", 
        "lasso2d",
        "zoomIn2d", 
        "zoomOut2d", 
        "sendDataToCloud",
        "toggleSpikelines" ,
        "logo"
        ],
    'displaylogo': False,
    }

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
    plot_colors =  [
        "#DDDDDD","#D6D9DC","#CFD6DB","#C8D3DA","#C1D0D9",
        "#BACDD9","#B3C9D8","#ACC6D7","#A5C3D6","#9EC0D6",
        "#97BDD5","#90BAD4","#89B6D3","#82B3D3","#7BB0D2",
        "#74ADD1","#6DAAD0","#66A6CF","#5FA3CF","#58A0CE",
        "#539DCC","#4F99CA","#4C95C8","#4992C6","#458EC3",
        "#428AC1","#3F87BF","#3B83BD","#3880BA","#347CB8",
        "#3178B6","#2E75B4","#2A71B1","#276DAF","#236AAD",
        "#2066AB","#1D62A8","#195FA6","#165BA4","#1358A2"
        ]
    levels = len(plot_colors)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return plot_colors[ind]

def plot_jaccard_knee_frag(count_data_file,outdir):
    count_data = pd.read_csv(count_data_file, index_col=0, sep=',')
    count_data['New']=count_data.index
    cell_bc = np.array(count_data[count_data['Beads'] == 'true'].index)
    sorted_bc = np.array(count_data.index)
    sorted_counts = np.array(count_data['UMI'])
    total_bc = len(sorted_bc)
    ix1 = count_data.drop_duplicates('Beads',keep='first').index[1]-1
    ix2 = count_data.drop_duplicates('Beads',keep='last').index[0]
    plot_segments = []
    barcodeSegment = collections.namedtuple(
        'barcodeSegment', 
        ['start', 'end', 'density', 'legend']
        )
    plot_segments.append(
        barcodeSegment(
            start=0, 
            end=ix1, 
            density=1.0, 
            legend=True
            )
        )
    plot_segments.append(
        barcodeSegment(
            start=ix2+1, 
            end=total_bc, 
            density=0.0, 
            legend=True
            )
        )
    mixed_segments = segment_log_plot(sorted_counts, ix1, ix2)
    for i in range(len(mixed_segments) - 1):
        num_cells = sum([1 for i in range(mixed_segments[i], mixed_segments[i + 1]) if sorted_bc[i] in cell_bc])
        density = float(num_cells)/float(mixed_segments[i + 1]-mixed_segments[i])
        plot_segments.append(
            barcodeSegment(
                start=mixed_segments[i], 
                end = mixed_segments[i + 1], 
                density=density, 
                legend=False
                )
            )
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
        
        data_dict = {
            "x": x,"y": y,"name": name,
            "hoverinfo": "text",
            "text": hover,
            "type": "scattergl",
            "mode": "lines",
            "line": {
                "width": 3,
                "color": plot_cmap(plot_segment.density),
                },
            "showlegend": plot_segment.legend,
            }
        plot_data.append(data_dict)

    plotly_data = [
        go.Scatter(
            x=dat['x'], 
            y=dat['y'], 
            name=dat['name'], 
            mode=dat['mode'], 
            showlegend=dat['showlegend'],
            marker={
                'color': dat['line']['color']
                }, 
            line=dat['line'], 
            text=dat['text']
            ) for dat in plot_data
            ]
    layout = go.Layout(
        xaxis = dict(
            type="log", 
            gridcolor="lightgrey",
            title="Barcode in Rank-descending Order",
            color="black",
            showline=True,
            zeroline=True,
            linewidth=1,
            fixedrange= True,
            linecolor="black"
            ),
        yaxis = dict(
            type="log",
            title="UMI counts",
            gridcolor="lightgrey",
            linewidth=1,
            fixedrange= True,
            color="black",
            linecolor="black"
            ),
        height=360,
        width=450,
        plot_bgcolor='rgba(0,0,0,0)',
        hovermode='closest',
        paper_bgcolor='white',
        legend = dict(
            x=1,
            y=1,
            traceorder="normal",
            font = dict(
                family="Arial",
                size=12,
                color="black"
                ),
            bordercolor="Black",
            borderwidth=0
            ),
        margin = dict(
            l=0,r=0,b=0,t=0,pad=1
            ),
        font = dict(
            size=10
            )
        )
    fig = go.Figure(
        data=plotly_data, 
        layout=layout
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/barcode_rank.html",
        auto_open=False,
        config=config
        )
    fig2 = py.offline.plot(
        fig, 
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir + "/barcode_rank.div",'w')
    fw.write(fig2)
    fw.close()

#### cluster
def plot_cluster(cluster_file,outdir):
    cluster = pd.read_csv(cluster_file)
    cluster[['Cluster']] = cluster[['Cluster']].astype('str')
    fig = px.scatter(
        cluster, 
        x=cluster.UMAP_1, 
        y=cluster.UMAP_2, 
        color= cluster['Cluster'],
        color_discrete_sequence=px.colors.qualitative.G10
        )
    
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
    
    fig.update_traces(
        marker={'size': 3}
        )
    fig.update_xaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/cluster.html",
        auto_open=False,
        config=config
        )
    fig2=py.offline.plot(
        fig, 
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir+"/cluster.div",'w')
    fw.write(fig2)
    fw.close()
    barplot_chsize = open(outdir+"/cluster.div","r").read()
    barplot_chsize=barplot_chsize.replace('width:100%','width:565px')
    fw1 = open(outdir+"/cluster_chsize.div","w")
    fw1.write(barplot_chsize)
    fw1.close()

### umi
def plot_cluster_umi(cluster_file,outdir):
    cluster = pd.read_csv(cluster_file)
    fig = px.scatter(
        cluster, 
        x=cluster.UMAP_1, 
        y=cluster.UMAP_2, 
        color= cluster['nUMI'],
        )
    
    fig.update_layout(
        autosize=False,
        width=520,
        height=500,
        legend_title=dict(font=dict(size=13),text='nUMI',),
        legend=dict(font=dict(size=10,family='Arial'),),
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',),
        yaxis=dict(gridcolor='lightgray',)
        )
    
    fig.update_traces(
        marker={'size': 3}
        )
    fig.update_xaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/cluster_nUMI.html",
        auto_open=False,
        config=config
        )
    fig2=py.offline.plot(
        fig, 
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir+"/nUMI.div",'w')
    fw.write(fig2)
    fw.close()
    barplot_chsize = open(outdir+"/nUMI.div","r").read()
    barplot_chsize=barplot_chsize.replace('width:100%','width:480px')
    fw1 = open(outdir+"/nUMI_chsize.div","w")
    fw1.write(barplot_chsize)
    fw1.close()

## cluster annotation plot
def plot_clusteranno(cluster_file,outdir):
    cluster = pd.read_csv(cluster_file)
    cluster = cluster.sort_values(by=['Predict_number','Cluster'],ascending=False)
    fig = px.scatter(
        cluster, 
        x=cluster.UMAP_1, 
        y=cluster.UMAP_2, 
        color= cluster['Predicted cell type'],
        color_discrete_sequence=px.colors.qualitative.G10
        )
    
    fig.update_layout(
        autosize=False,
        width=800,
        height=500,
        legend_title=dict(font=dict(size=16),text='Predicted cell type: cell number',),
        legend=dict(x=1.2,y=0.5,font=dict(size=10,family='Arial'),itemsizing='constant'),
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',),
        yaxis=dict(gridcolor='lightgray',)
    )

    fig.update_traces(
        marker={'size': 3}
        )
    fig.update_xaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/anno.html",
        auto_open=False,
        config=config
        )
    fig2=py.offline.plot(
        fig, 
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir+"/anno.div",'w')
    fw.write(fig2)
    fw.close()
    barplot_chsize = open(outdir+"/anno.div","r").read()
    barplot_chsize=barplot_chsize.replace('width:100%','width:565px')
    fw1 = open(outdir+"/anno_chsize.div","w")
    fw1.write(barplot_chsize)
    fw1.close()


def plot_saturation(saturation_file,outdir):
    saturation = pd.read_csv(saturation_file, sep="\t",header=0)
    x=saturation['Mean Reads per Cell']
    y=saturation['Sequencing Saturation']
    if len(saturation) > 2:
        xnew = np.linspace(x.min(),x.max(),50)
        ynew = make_interp_spline(x,y)(xnew)
    else:
        xnew = x
        ynew = y
    fig = px.line(saturation, x=xnew, y=ynew )
    fig.update_layout(
        autosize=False,
        width=565,
        height=500,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',title="Mean Reads per Cell"),
        yaxis=dict(gridcolor='lightgray',title="Sequencing Saturation"),
        yaxis_range=[0,100]
    )
    fig.update_xaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_traces(
        line=dict(color="#337ab7", width=3)
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/saturation.html",
        auto_open=False,
        config=config
        )
    fig2 = py.offline.plot(
        fig, 
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir+"/saturation.div",'w')
    fw.write(fig2)
    fw.close()


    x=saturation['Mean Reads per Cell']
    y=saturation['Median Genes per Cell']
    if len(saturation) > 2:
        xnew = np.linspace(x.min(),x.max(),50)
        ynew = make_interp_spline(x,y)(xnew)
    else:
        xnew = x
        ynew = y
    fig = px.line(saturation, x=xnew, y=ynew )
    fig.update_layout(
        autosize=False,
        width=520,
        height=500,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',title="Mean Reads per Cell"),
        yaxis=dict(gridcolor='lightgray',title="Median Genes per Cell")
    )
    fig.update_xaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_traces(
        line=dict(color="#337ab7", width=3)
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/saturation2.html",
        auto_open=False,
        config=config
        )
    fig2 = py.offline.plot(
        fig, 
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir+"/saturation2.div",'w')
    fw.write(fig2)
    fw.close()


### plot log_uniqueFrags
def plot_cluster_uniqueFrags(cluster_file,outdir):
    cluster = pd.read_csv(cluster_file)
    fig = px.scatter(
        cluster, 
        x=cluster.UMAP_1, 
        y=cluster.UMAP_2, 
        color= cluster['log10_uniqueFrags'],
        color_continuous_scale=px.colors.sequential.Viridis
        )
    
    fig.update_layout(
        autosize=False,
        width=520,
        height=500,
        legend_title=dict(font=dict(size=13),text='log10_uniqueFrags',),
        legend=dict(font=dict(size=10,family='Arial'),),
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',),
        yaxis=dict(gridcolor='lightgray',)
        )
    
    fig.update_traces(
        marker={'size': 3}
        )
    fig.update_xaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True, 
        zerolinewidth=1, 
        zerolinecolor='gray'
        )
    py.offline.plot(
        fig, 
        filename = outdir+"/cluster_lg_uniqueFrags.html",
        auto_open=False,
        config=config
        )
    fig2=py.offline.plot(
        fig,
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )
    fw = open(outdir+"/lg_uniqueFrags.div",'w')
    fw.write(fig2)
    fw.close()
    barplot_chsize = open(outdir+"/lg_uniqueFrags.div","r").read()
    barplot_chsize=barplot_chsize.replace('width:100%','width:480px')
    fw1 = open(outdir+"/lg_uniqueFrags_chsize.div","w")
    fw1.write(barplot_chsize)
    fw1.close()

### atac
def plot_barcode_atac_frag(barcodeCount,d2cCutoff,outdir):
    seqnum = pd.read_csv(barcodeCount,encoding="utf-8",header=None,sep="\t")
    seqnum.columns = ['barcode','num']
    seqnum_sort = seqnum.sort_values(by='num',ascending=False)
    # seqnum_sort = seqnum_sort[seqnum_sort['num']>=50]
    # new_data = pd.DataFrame({'num': list(range(1, 11)),'barcode': ['barcode']*10})
    # seqnum_sort = pd.concat([seqnum_sort, new_data], ignore_index=True)
    seqnum_sort['range'] = range(1,len(seqnum_sort)+1)
    seqnum_sort_first = set(seqnum_sort[seqnum_sort[["num"]].duplicated(keep="first")].index)
    seqnum_sort_last = set(seqnum_sort[seqnum_sort[["num"]].duplicated(keep="last")].index)
    seqnum_sort_inter = seqnum_sort_first & seqnum_sort_last
    seqnum_sort=seqnum_sort.drop(list(seqnum_sort_inter ),axis=0)

    df2 = pd.read_csv(d2cCutoff,encoding="utf-8",header=None,sep="\t")
    beads_number = int(df2[1][0])
    seqnum_sort_True = seqnum_sort[(seqnum_sort['num'] >= beads_number)]
    seqnum_sort_False = seqnum_sort[(seqnum_sort['num'] < beads_number)]
    trace0_x = list(seqnum_sort_True['range'])
    trace0_y = list(seqnum_sort_True['num'])  
    trace1_x = list(seqnum_sort_False['range'])
    trace1_y = list(seqnum_sort_False['num'])
    blue_line = list(zip(trace0_x, trace0_y))
    blue_line = [list(i) for i in blue_line]
    red_line = list(zip(trace1_x, trace1_y))
    red_line = [list(i) for i in red_line]

    trace0 = go.Scattergl(
        x = trace0_x,
        y = trace0_y,
        mode="lines",
        name="TRUE",
        line=dict(color="#005bac",width=5)
    )
    trace1 = go.Scattergl(
        x = trace1_x,
        y = trace1_y,
        mode="lines",
        name="FALSE",
        line=dict(color="gray",width=5)
    )
    config={'displayModeBar': False}
    #hovermode='closest',)
    data = [trace0, trace1]
    config={'displayModeBar': False}
    layout = go.Layout(
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
                        title="Unique Fragments",
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
    data = [trace0, trace1]
    fig = dict(data=data, layout=layout)
    py.offline.plot(fig, filename = outdir+"/barcode_rank.html",auto_open=False,config=config)
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(outdir+"/barcode_rank.div",'w')
    fw.write(fig2)

def plot_jaccard_atac_frag(CorrelationBarcodes,d2cCutoff,outdir):
    df = pd.read_csv(CorrelationBarcodes,encoding="utf-8",compression="gzip",sep="\t") 
    jaccard = np.array(df['jaccard_distance']).tolist()
    jaccard.sort(key=float) 
    jaccard = [i for i in jaccard if i > 0.001]
    ls_len = len(jaccard)
    sort_num = [i for i in range(ls_len+1)]
    sort_num.sort(reverse=True) 
    cell_info = pd.read_csv(d2cCutoff,encoding="utf-8",header=None,sep="\t")
    trans_point = round(float(cell_info[1][1]),5)
    jaccard_sort = list(zip(sort_num,jaccard))
    higher_than_point_trans = [i for i in jaccard_sort if float(i[1]) > trans_point]
    lower_than_point_trans = [i for i in jaccard_sort if float(i[1]) <= trans_point]
    trace0_x = [i[0] for i in higher_than_point_trans]
    trace0_x.sort()  
    trace0_y = [i[1] for i in higher_than_point_trans]
    trace0_y.sort(reverse=True)
    trace1_x = [i[0] for i in lower_than_point_trans]
    trace1_x.sort()
    trace1_y = [i[1] for i in lower_than_point_trans]
    trace1_y.sort(reverse=True) 
    blue_line = list(zip(trace0_x, trace0_y))
    blue_line = [list(i) for i in blue_line]
    black_line = list(zip(trace1_x, trace1_y))
    black_line = [list(i) for i in black_line] 
    trace0 = go.Scatter(
        x = trace0_x,
        y = trace0_y,
        mode="lines",
        name="TRUE",
        line=dict(color="#005bac",width=5)
    )
    trace1 = go.Scatter(
        x = trace1_x,
        y = trace1_y,
        mode="lines",
        name="FALSE",
        line=dict(color="grey",width=5)
        
    )
    config={'displayModeBar': False}
    layout = go.Layout(
                        xaxis=dict(type="log", 
                        gridcolor="lightgrey",
                        title="D2C Overlap Score in Rank-descending Order",
                        color="black",
                        showline=True,
                        zeroline=True,
                        linewidth=1,fixedrange= True,
                        linecolor="black"
                        ),
                        
                        yaxis = dict(
                        type="log",
                        title="D2C Overlap Score per Barcode Pair",
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
                        bordercolor="black",
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
    data = [trace0, trace1]
    fig = dict(data=data, layout=layout)
    py.offline.plot(fig, filename = outdir +"/jaccard_rank.html",auto_open=False,config=config) 
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config) 
    fw = open(outdir +"/jaccard_rank.div",'w') 
    fw.write(fig2)

def atac_saturation(sequenceSaturation,outdir):
    saturation=pd.read_csv(sequenceSaturation,sep="\t")
    x=saturation['mean_frags_per_cell']
    y=saturation['saturation']
    fig = px.line(saturation, x=x, y=y)

    fig.update_layout(
        autosize=False,
        width=500,
        height=400,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',title="Mean Fragments per Cell"),
        yaxis=dict(gridcolor='lightgray',title="Sequencing Saturation"),
    )
    fig.update_xaxes(
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='gray'
        )
    fig.update_yaxes(
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='gray'
        )
    fig.update_traces(
        line=dict(color="#337ab7", width=3)
        )
    py.offline.plot(
        fig,
        filename = outdir+"/saturation.html",
        auto_open=False,
        config=config
        )
    fig2 = py.offline.plot(
        fig,
        include_plotlyjs=False,
        show_link=False,
        output_type='div',
        config=config
        )

    # data=go.Scattergl(x=df["mean_frags_per_cell"],y=df["saturation"], mode='lines',line=dict(color="#337ab7", width=3),)
    # layout = go.Layout(xaxis=dict(
    #                     gridcolor="lightgrey",
    #                     title="mean fragments per cell",
    #                     color="black",
    #                     showline=True,
    #                     zeroline=True,
    #                     linewidth=1,fixedrange= True,
    #                     linecolor="black"
    #                     ),
    #                 yaxis = dict(
    #                     title="saturation",
    #                     gridcolor="lightgrey",
    #                     linewidth=1,fixedrange= True,
    #                     color="black",
    #                     linecolor="black"
    #                     ),
    #                 height=450,width=500,
    #                 plot_bgcolor='rgba(0,0,0,0)',
    #                 hovermode='closest',
    #                 paper_bgcolor='white',
    #                 #title="saturation : "+str(round(float(df['saturation'].tolist()[-1]),3))
    # )
    # fig=dict(data=data, layout=layout)
    # config={'displayModeBar': False}
    # py.offline.plot(fig, filename =outdir+"/saturation.html")
    # fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(outdir+"/saturation.div",'w')
    fw.write(fig2)
    fw.close()