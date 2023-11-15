# -*- coding: utf-8 -*-
"""
Created on Fri June 18 14:33:50 2021

@author: qwzhou
=======================================
plot boxplot and etc
=======================================
This is a script for basic plot of batmeth2 output
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import pandas as pd 
import matplotlib.colors as pltc
all_colors = [k for k,v in pltc.cnames.items()]
import seaborn as sns

#plt.switch_backend('agg')

def writableFile(string):
    """
    function that tests if a given path is writable
    """
    try:
        open(string, 'w').close()
        os.remove(string)
    except:
        msg = "{} file can't be opened for writing".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string

def getArgs(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-f", "--mrfile", 
                        default='', 
                        help="DNA AverMethylevel files, seperate by space. eg. wildtype.body/bin.txt", 
                        nargs='+')
    parser.add_argument("-l", "--label", default='wildtype', 
                        help="Labels of samples, sperate by space. eg. -l CG CHG", 
                        nargs='+')
    parser.add_argument("--chr_select", default='all', 
                        help="chromosome selected for boxplot and hexsin, default all. all means used all.", 
                        )
    parser.add_argument("-p", "--plotmode", default='boxnplot',
                        help="plot mode, [boxplot, boxnplot, all] default boxnplot.",
                        )
    parser.add_argument("--plotpoint", default='N',
                        help="plot point inside boxplot or boxnplot, defalt N.",
                        )
    parser.add_argument("--pointsize", default=1,
                        choices=[1,2,3,4],
                        help="point size while plot point.",
                        type=int)

    ## coverfile
    parser.add_argument("-c", "--coverfile", 
                        default='', 
                        help="DNA AverMethylevel files, seperate by space. eg. wildtype.NCcoverage.txt", 
                        nargs='+')
    # me stats
    parser.add_argument("-s", "--statsfile",
                        default='',
                        help="DNA AverMethylevel stats file",
                        nargs='+')
    parser.add_argument('--outfilename', '-o',
                       help='Output file name prefix.',
                       metavar='FILENAME',
                       default='myMeth',
                       type=writableFile,
                       required=True)
    parser.add_argument('--yMin',
                        default=[0.0],
                        nargs='+',
                        help='Minimum value for the Y-axis. Multiple values, separated by '
                            'spaces can be set for each profile. If the number of yMin values is smaller than'
                            'the number of plots, the values are recycled.',
                        type=float)
    parser.add_argument('--yMax',
                        default=[],
                        nargs='+',
                        help='Maximum value for the Y-axis. Multiple values, separated by '
                            'spaces can be set for each profile. If the number of yMin values is smaller than'
                            'the number of plots, the values are recycled.',
                        type=float)
    parser.add_argument('--legend',
                        default=0,
                        choices=[0,1,2,3,4,5,6,7,8,9,10,11,12],
                        help='The location of the legend. '
                        'best	     :  0, '
                        'upper right :	1, '
                        'upper left  :  2, '
                        'lower left  :	3, '
                        'lower right :	4, '
                        'right	     :  5, '
                        'center left :  6, '
                        'center right:	7, '
                        'lower center:	8, '
                        'upper center:	9, '
                        'center	     : 10, '
                        'out         : 11, '
                        'none        : 12',
                        type=int)
    parser.add_argument('--legendsize',
                        default=10,
                        help='the text size of the legend.',
                        type=int
                        )
    parser.add_argument("-ft", "--image_format",
                        default="",
                        help="The file format, e.g. 'png', 'pdf', 'svg', ... "
                        "The behavior when this is unset is documented under fname.")
    parser.add_argument("--dpi",
                        default=200,
                        help="Set the DPI to save the figure. default: 200",
                        type=int)
    parser.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
                        
    return parser

import collections

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

#a549.bt2.NCcoverage.txt
#[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
def readcoverfile(coverfile, xdata, NCcovergae):
    with open(coverfile, 'r') as mdata:
        nline=0
        for line in mdata:
            line=line.strip('\n')
            Ncover=line.split()
            nline+=1
            if nline == 1 and len(xdata)==0 :
                xdata.append(Ncover[1:])
            elif nline>1:
                NCcovergae.append(Ncover[1:])
    #print(NCcovergae[0])

def readstatsfile(cfile, Sxdata, Sydaya, Slydata):
    with open(cfile, 'r') as mdata:
        nline=0
        for line in mdata:
            line=line.strip('\n')
            Nstats=line.split()
            nline+=1
            if nline == 1 and len(Sxdata)==0 :
                Sxdata.append(Nstats[1:])
            elif nline>1:
                #if len(Sydaya) == 0:
                #    for i in range(5):
                #        j=[]
                #        Sydaya.append(j)
                #for i in range(0, len(Nstats), 1):
                #    Sydaya[i].append(Nstats[i])
                Sydaya.append(list(map(eval, Nstats[1:])))
                Slydata.append(list(map(eval, Nstats[1:])))

def plotbarh(results, category_names, outfilename, mydpi, image_format, legend, legendsize):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.colormaps['RdYlGn_r'](
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(5.2, 3))
    ax.invert_yaxis()
    #ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    print('size', len(category_names), len(category_colors))
    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        print((colname,color))
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        print("widths", widths, "ccc", data_cum)

        print('xxx', labels, widths, starts)
        rects = ax.barh(labels, widths, left=starts, height=0.5,
                        label=colname, color=color)
        #rects = ax.bar(labels, widths, bottom=starts, width=0.5,
        #                label=colname, color=color)

        # number label
        #r, g, b, _ = color
        #text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        #ax.bar_label(rects, label_type='center', color=text_color)

    #if legend < 11:
    #    plt.legend(loc=legend, prop={'size': legendsize}, frameon=False)
    #elif legend == 11:
    #    plt.legend(bbox_to_anchor=(1.01, 0.05), loc=3, borderaxespad=0, prop={'size': legendsize}, frameon=False)

    plt.xlabel("Count")
    plt.ylabel("Samples")
    plt.title("DNA methylation level category")

    plt.legend(ncol=5, bbox_to_anchor=(0.98, -0.02), loc=4, frameon=False, bbox_transform=fig.transFigure,
        fontsize='small') #ncol=len(category_names), bbox_to_anchor=(0, 1),

    fig.tight_layout()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)
    #return fig, ax


def plotbar(xlabels, ydata, lydata, samplelabels, outfilename, mydpi, image_format, legend, legendsize):
    #plt.style.use('fivethirtyeight')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    # 此处的 _ 下划线表示将循环取到的值放弃，只得到[0,1,2,3,4]
    xlabel = [xlabels[0]]
    ind = [x for x, _ in enumerate(xlabel)] 
    # line plot
    #plt.plot(ind, y_data, linestyle="--")
    #plt.plot(ind, y2_data+y_data, linestyle="-.")
    Ny = len(ydata[0])
    print(ind, Ny, xlabel)
    #if Ny > 1:
    #    for i in range(0,Ny,1):
    #        plt.plot(ind, lydata[i], linestyle="--")
    Ny = len(ydata)
    #绘制堆叠图
    #plt.bar(ind, bronzes, width=0.5, label='bronzes', width=0.67)
    #plt.bar(ind, silvers, width=0.5, label='silvers', color='silver', width=0.67, bottom=bronzes)
    #plt.bar(ind, golds, width=0.5, label='golds', color='gold', width=0.67, bottom=silvers+bronzes)
    #mybottom = ydata[0]
    #for i in range(0,Ny,1):
    #    if i == 0:
    #        plt.bar(ind, ydata[i], width=0.6, label=samplelabels[i])
    #    else:
    #        for j in range(1,i-1,1):
    #            mybottom = mybottom + ydata[i]
    #        plt.bar(ind, ydata[i], width=0.6, label=samplelabels[i], bottom=mybottom)
    print('ydaya', ydata)
    df = pd.DataFrame(ydata, index=samplelabels, columns=['a', 'b', 'c', 'd', 'e'])
    print(df.index, samplelabels, xlabel, ydata[0], "sam", df, df['a'], "EEE")
    ax.bar(df.index, df['a'], width=0.6, label="[0-0.2)")
    ax.bar(df.index, df['b'], width=0.6, label="[0.2-0.4)", bottom=df['a'])
    ax.bar(df.index, df['c'], width=0.6, label="[0.4-0.6)", bottom=df['a']+df['b'])
    ax.bar(df.index, df['d'], width=0.6, label="[0.6-0.8)", bottom=df['a']+df['b']+df['c'])
    ax.bar(df.index, df['e'], width=0.6, label="[0.8-1)", bottom=df['a']+df['b']+df['c']+df['d'])
    #设置坐标轴
    if legend < 11:
        plt.legend(loc=legend, prop={'size': legendsize}, frameon=False)
    elif legend == 11:
        plt.legend(bbox_to_anchor=(1.01, 0.05), loc=3, borderaxespad=0, prop={'size': legendsize}, frameon=False)
    #plt.xticks(df.index, xlabel) 
    plt.ylabel("Count") 
    plt.xlabel("Samples") 
    plt.legend(loc="upper right") 
    plt.title("DNA methylation level category")
    fig.tight_layout()
    #plt.show()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)

'''
line style
'-' or 'solid'	solid line
'--' or 'dashed'	dashed line
'-.' or 'dashdot'	dash-dotted line
':' or 'dotted'	dotted line
'None' or ' ' or ''	draw nothing
'''

'''
marker
'.'	point marker
','	pixel marker
'o'	circle marker
'v'	triangle_down marker
'^'	triangle_up marker
'<'	triangle_left marker
'>'	triangle_right marker
'1'	tri_down marker
'2'	tri_up marker
'3'	tri_left marker
'4'	tri_right marker
'8'	octagon marker
's'	square marker
'p'	pentagon marker
'P'	plus (filled) marker
'*'	star marker
'h'	hexagon1 marker
'H'	hexagon2 marker
'+'	plus marker
'x'	x marker
'X'	x (filled) marker
'D'	diamond marker
'd'	thin_diamond marker
'|'	vline marker
'_'	hline marker
'''
import matplotlib.ticker as ticker
legendsize=10
def plotline11(x, NCcovergae, outfilename, mydpi, image_format, legend, label):
    x = list(map(eval, x[0]))
    print('xxx', x)
    figextend = 0
    if legend == 11:
        figextend = 2
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5+figextend, 4))

    colors= ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    linstyle = ['-', '--', '-.', ':']
    marker = ['o', 'v', '^', 'x', '*', '+', '1', '2', '>', '<', '3', '4', 'h', 'H']
    for idx in range(0, len(NCcovergae), 1):
        y = NCcovergae[idx] #list(map(eval, NCcovergae[idx]))
        ax.plot(x, y, color=colors[idx] ,marker='o', linestyle='dashed', label=label[idx],
          linewidth=1, markersize=3)
    #x_major_locator=ticker.MultipleLocator(1)
    #ax.xaxis.set_major_locator(x_major_locator)
    ax.set_title("The number of C in bins", fontsize=16)
    ax.set_xlabel("Methy level", fontsize=14)
    ax.set_ylabel("Number of cytocines", fontsize=14)

    if legend < 11:
        plt.legend(loc=legend, prop={'size': legendsize}, frameon=False, labels=label)
    elif legend == 11:
        plt.legend(bbox_to_anchor=(1.01, 0.05), loc=3, borderaxespad=0, prop={'size': legendsize}, frameon=False, labels=label)

    fig.tight_layout()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)

def plotline(x, NCcovergae, outfilename, mydpi, image_format, legend, label):
    #x = np.linspace(0, 10, 500)
    #y = np.sin(x)
    #x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    #y = list(map(eval, NCcovergae[0]))
    #fig, ax = plt.subplots()
    x = list(map(eval, x[0]))
    print('xxx', x)
    figextend = 0
    if legend == 11:
        figextend = 2
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9+figextend, 4))

    colors= ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    linstyle = ['-', '--', '-.', ':']
    marker = ['o', 'v', '^', 'x', '*', '+', '1', '2', '>', '<', '3', '4', 'h', 'H']
    for idx in range(0, len(NCcovergae)-1, 2):
        y = list(map(eval, NCcovergae[idx]))
        ax[0].plot(x, y, color=colors[int(idx/2)] ,marker='o', linestyle='dashed', label=label[int(idx/2)],
          linewidth=1, markersize=3)
    x_major_locator=ticker.MultipleLocator(1)
    ax[0].xaxis.set_major_locator(x_major_locator)
    ax[0].set_title("Coverage distribution of cytosines", fontsize=16)
    ax[0].set_xlabel("Coverage threshold (>=)", fontsize=14)
    ax[0].set_ylabel("Number of cytocines", fontsize=14)

    # Using plot(..., dashes=...) to set the dashing when creating a line
    ymin=1
    for idx in range(1, len(NCcovergae), 2):
        y = list(map(eval, NCcovergae[idx]))
        ax[1].plot(x, y, marker=marker[int(idx/2)], linestyle='dashed', linewidth=1) #, label='Coverage percent' 
        ymin = min(y) if min(y)<ymin else ymin
    ymin = ymin - 0.02
    if ymin < 0:
        ymin = 0
    x_major_locator=ticker.MultipleLocator(1)
    ax[1].xaxis.set_major_locator(x_major_locator)
    ax[1].set_title("Coverage distribution of cytosines", fontsize=16)
    ax[1].set_xlabel("Coverage threshold (>=)", fontsize=14)
    ax[1].set_ylabel("Percent of cytosines", fontsize=14)
    ax[1].set(ylim=(ymin, 1.005))
    if legend < 11:
        plt.legend(loc=legend, prop={'size': legendsize}, frameon=False, labels=label)
    elif legend == 11:
        plt.legend(bbox_to_anchor=(1.01, 0.05), loc=3, borderaxespad=0, prop={'size': legendsize}, frameon=False, labels=label)
    
    #ax[0].legend()
    #ax[1].legend()
    fig.tight_layout()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)
    #plt.show()

dict_data = {}
def readmethfile(methlevel, methleveldf, filename, label, chr_choice, nsample, processed):
    # chr pos strand mC mCT ID
    global dict_data
    with open(filename, 'r') as mdata:
        for line in mdata:
            data = line.split()
            chrom, pos, strand, mC, cover, lineID = data
            if chr_choice != 'all' and chrom != chr_choice:
                continue
            mr = int(mC)/int(cover)
            methlevel.setdefault(label,[]).append(float(mr))
            ## df
            methleveldf.setdefault("meth",[]).append(float(mr))
            methleveldf.setdefault("context",[]).append(label)
            dict_data.setdefault(chrom+':'+pos,[]).append(float(mr))

    if(processed+1<nsample):
        return
    if(nsample>1):
        #for key, val in dict_data.items():
        for key in list(dict_data.keys()):
            if(len(dict_data[key])<nsample):
                del dict_data[key]

    if(len(dict_data)==0):
        print("The merged matrix data is zero")
        exit()
    frame = pd.DataFrame.from_dict(dict_data,orient = 'index')
    return frame

'''
sns.set
context='' 参数控制着默认的画幅大小，分别有 {paper, notebook, talk, poster} 四个值。其中，poster > talk > notebook > paper。
style='' 参数控制默认样式，分别有 {darkgrid, whitegrid, dark, white, ticks}，你可以自行更改查看它们之间的不同。
palette='' 参数为预设的调色板。分别有 {deep, muted, bright, pastel, dark, colorblind} 等，你可以自行更改查看它们之间的不同。
剩下的 font='' 用于设置字体，font_scale= 设置字体大小，color_codes= 不使用调色板而采用先前的 'r' 等色彩缩写。
'''
#import seaborn as sns
def plotboxsns(mdata, xl, yl,orient):
    plt.figure(figsize=(6,5))        #定义图像大小
    mdatadf = pd.DataFrame.from_dict(mdata)
    #sns.set(context='paper',font_scale=2,style='white',
    #       palette='colorblind')
    #plt.rc('font',family='Times New Roman',size=12)
    if orient == 'v':
        labels = mdatadf[xl].unique()
        Nlabel = len(labels)
    else:
        labels = mdatadf[yl].unique()
        Nlabel = len(labels)
    ax = sns.boxplot(data=mdatadf, x=xl, y=yl,        #传入数据
                linewidth=3,        #箱边线宽度
                width=0.8,        #箱体宽度，默认0.8
                whis=1.5,        #计算上限和下限时四分位距(IQR)前的系数，默认1.5
                showfliers=True,        #是否显示异常值
                fliersize=5,        #异常值大小，默认5
                #hue="context", # legend
                palette="Set3",        #颜色盘
                #order=class_order,        #顺序
                orient=orient, # v h
                notch=False,
                saturation=1)        #颜色饱和度，默认0.75
    
    g2 = sns.swarmplot(data=mdatadf, x=xl, y=yl,
              color='grey', size=4, linewidth=0.5, edgecolor='k',
              #order=class_order, 
              alpha=0.75)
    
    for i in range(len(ax.artists)):
        box = ax.artists[i]
        box.set_edgecolor('red')
        box.set_facecolor('white')
        for j in range(5):
            k = i*5 + j
            line = ax.lines[k]
            line.set_color('red')

    '''
      Horizontal alignment must be one of ('center', 'right', 'left')
	  Vertical alignment must be one of ('top', 'bottom', 'center', 'baseline', 'center_baseline')
    '''
    #plt.xlim([0,12])  # x轴边界
    #plt.ylim([0,1.5])  # y轴边界
    if orient == 'v':
        ylabel = 'DNA methylation'
        plt.ylabel(ylabel, fontsize=16)
        plt.yticks(fontsize=12)
        plt.xlabel('')
        #ax.set_xticklabels(ax.get_xticklabels(), 
        #         rotation=60,
        #         ha="center",
        #         va='bottom',
        #         fontsize=15,
        #         rotation_mode='anchor'
        #         )
        plt.xticks(ticks=range(Nlabel), labels=labels,
                rotation=60,        #设置刻度文字旋转角度
                fontsize=14,
                ha='right', va='center',        #刻度文字对齐方式，当rotation_mode为’anchor'时，对齐方式决定了文字旋转的中心。ha也可以写成horizontalalignment，va也可以写成verticalalignment。
                rotation_mode='anchor')
    else:
        xlabel = 'DNA methylation'
        plt.xlabel(xlabel, fontsize=16)
        plt.xticks(fontsize=12)
        plt.ylabel('')
        #ax.set_yticklabels(ax.get_yticklabels(), 
        #         rotation=60,
        #         ha="center",
        #         va='bottom',
        #         fontsize=15,
        #         rotation_mode='anchor'
        #         )
        plt.yticks(ticks=range(Nlabel), labels=labels,
                fontsize=14, ha='center', va='bottom',
                rotation=60, 
                rotation_mode='anchor')

    
    #ax = plt.axes()
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)

    plt.show()
    #plt.savefig(outfig,
    #        dpi=300,        # 设置分辨率
    #        format=None,        #设置图片格式，默认None，如果未设置使用文件名设置的格式 eps, pdf, pgf, png, ps, raw, rgba, svg, svgz。
    #        bbox_inches='tight',        #设置为tight，防止有时图片保存不完整
    #        facecolor='w',        #背景颜色，默认'w'白色
    #        edgecolor='w')        #边框颜色，默认'w'白色

def plotviolionsns(mdata, xl, yl,orient):
    plt.figure(figsize=(6,5))        #定义图像大小
    mdatadf = pd.DataFrame.from_dict(mdata)

    if orient == 'v':
        labels = mdatadf[xl].unique()
        Nlabel = len(labels)
    else:
        labels = mdatadf[yl].unique()
        Nlabel = len(labels)
    #sns.set(font_scale = 2) #字体大小
    g = sns.violinplot(data=mdatadf,         #传入数据，对sepal_length这一列画图，根据class分组
                   x=xl,
                   y=yl,
                   linewidth=3, # 线宽
                   inner='box',  # 内部数据形式，默认为box，内部画个小箱型图 'box' quartiles None stick point
                   #order="",   # 指定顺序
                   #hue="context", # legend
                   bw=0.5, #
                   orient=orient, # v h
                   width=0.9, #width：float，控制钢琴图的宽度（比例）
                   scale='count', # scale：该参数用于缩放每把小提琴的宽度，有“area”, “count”, “width”三种方式
                   cut=2, # cut：float，距离，以带宽大小为单位，以控制小提琴图外壳延伸超过内部极端数据点的密度。设置为0以将小提琴范围限制在观察数据的范围内（即，在ggplot中具有与trim = true相同的效果)
                   saturation=1) 
    
    g2 = sns.swarmplot(data=mdatadf, x=xl, y=yl,
                   color='grey', alpha=0.5, size=4,        #颜色，透明度，大小
                   linewidth=0.5, edgecolor='black',        #边线宽度，边线颜色
                   #order=class_order, #顺序和violinplot保持一致
                   )   

    if orient == 'v':
        ylabel = 'Sepal length'
        plt.ylabel(ylabel, fontsize=18)        #设置y轴标签
        plt.yticks(fontsize=15)        #设置y轴刻度字体大小
        plt.xlabel('')        #去掉x轴标签
        plt.xticks(ticks=range(Nlabel),        #设置要显示的x轴刻度，若指定空列表则去掉x轴刻度
                labels=labels,        #设置x轴刻度显示的文字，要与ticks对应   
                fontsize=15,        #设置刻度字体大小
                rotation=60,        #设置刻度文字旋转角度
                ha='right', va='center',        #刻度文字对齐方式，当rotation_mode为’anchor'时，对齐方式决定了文字旋转的中心。ha也可以写成horizontalalignment，va也可以写成verticalalignment。
                rotation_mode='anchor')        #我的设置表示文字以右边线的中点为中心旋转。
    else:
        xlabel = 'Sepal length'
        plt.xlabel(xlabel, fontsize=18)        #设置y轴标签
        plt.xticks(fontsize=15)        #设置y轴刻度字体大小
        plt.ylabel('')        #去掉x轴标签
        plt.yticks(ticks=range(Nlabel),        #设置要显示的x轴刻度，若指定空列表则去掉x轴刻度
                labels=labels,        #设置x轴刻度显示的文字，要与ticks对应   
                fontsize=15,        #设置刻度字体大小
                rotation=60,        #设置刻度文字旋转角度
                ha='center', va='bottom',        #刻度文字对齐方式，当rotation_mode为’anchor'时，对齐方式决定了文字旋转的中心。ha也可以写成horizontalalignment，va也可以写成verticalalignment。
                rotation_mode='anchor')        #我的设置表示文字以右边线的中点为中心旋转。

    #ax = plt.axes()
    #ax.spines['top'].set_visible(False)        #去掉图像上边框和右边框
    #ax.spines['right'].set_visible(False)

    plt.show()

'''
point correlation
'''
def pointcorplot(mdata, outfilename, mydpi, image_format):
    from scipy.stats import gaussian_kde
    medatas = list(mdata.values())
    x = medatas[0]
    y = medatas[1]

    #cor
    method="spearman" #kendall
    xdf = pd.Series(x)
    ydf = pd.Series(y)
    mecor = round(xdf.corr(ydf, method),4)

    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = np.array(x)[idx], np.array(y)[idx], np.array(z)[idx]

    fig, ax = plt.subplots()
    plt.scatter(x, y, c=z, s=20, cmap='Spectral')
    ax.set_xlabel(method + ': ' + str(mecor), fontsize= 12)
    ax.set_ylabel("DNA methylation", fontsize= 12)
    ax.set_title("DNA methylation density")
    plt.colorbar()
    fig.tight_layout()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)
    #plt.show()

'''
correlation
'''
def corplot(mdata, outfilename, mydpi, image_format):
    medatas = list(mdata.values())
    x = medatas[0]
    y = medatas[1]
    xmin = min(x)
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)
    #methdf = pd.DataFrame(mdata)
    method="spearman" #kendall
    #print(round(methdf.corr(method),4))

    xdf = pd.Series(x)
    ydf = pd.Series(y)
    mecor = round(xdf.corr(ydf, method),4)

    fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(9, 4))
    fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
    ax = axs[0]
    hb = ax.hexbin(x, y, gridsize=50, cmap='viridis',
                  marginals = True,
                  linewidths = 1)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    #ax.text(0.2, 0.8, method + ': ' + str(mecor), rotation=45)
    ax.set_xlabel(method + ': ' + str(mecor), fontsize= 12)
    ax.set_title("DNA methylation hexbin")
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('DNA meth')

    ax = axs[1]
    hb = ax.hexbin(x, y, gridsize=50, bins='log', cmap='viridis',
                   marginals = True,
                   linewidths = 1)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    #ax.set_xlabel(method + ': ' + str(mecor), fontsize= 12)
    ax.set_title("DNA methylation hexbin log")
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('log10(me)')

    fig.tight_layout()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)
    #plt.show()

'''
box plot and violin plot
'''
#from beeswarm import *
def plotboxplot(fig, axs, mdata, mdatadf, notch, vert, showbox, showfliers, showmeans, showmedians, boxprops, medianprops,
              meanprops, capprops, whiskerprops, fontsize, plotpoint, pointsize):

    # generate some random test data
    #all_data = [np.random.normal(0, std, 100) for std in range(6, 10)]
    all_data = list(mdata.values())
    labels = list(mdata.keys())

    # plot box plot
    notch = notch
    pos = np.arange(len(labels))
    bplot = axs.boxplot(all_data,
                    patch_artist=True,  # fill with color
                    notch=notch,  # notch shape
                    vert=vert,  # vertical box alignment, True
                    showbox=showbox, # Show the central box, True
                    showfliers=showfliers, # Show the outliers beyond the caps, False
                    showmeans=showmeans, # Show the arithmetic means, False
                    boxprops=boxprops, # dict, The style of the box, {'color':'orangered','facecolor':'pink'}
                    whis=1.5, # float or (float, float), default: 1.5, The position of the whiskers.
                    medianprops=medianprops, # dict, The style of the median, {'color':'green','linewidth':'1.5'}
                    meanprops=meanprops, # dict, The style of the median
                    capprops=capprops, # dict, The style of the caps, {'color':'black','linewidth':'2.0'}
                    whiskerprops=whiskerprops, # dict, The style of the whiskers
                    sym = '', # symbol of outliers, * #
                    labels=labels,
                    positions=pos,
                    )
    axs.set_title('DNA methylation box plot', fontsize=14)
    # fill with colors
    #for bplot in (bplot1, bplot2):
    # add fill color
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('green') ## 边框
        patch.set_linestyle('--') #'-' or 'solid', '--' or 'dashed', '-.' or 'dashdot', ':' or 'dotted', ''
        patch.set_linewidth(1.5)
    for patch, color in zip(bplot['medians'], colors):
        patch.set_color('green')
        patch.set_linestyle('-.')
        patch.set_linewidth('1.5')
    for patch, color in zip(bplot['caps'], colors):
        patch.set_color('green')
        patch.set_linestyle('--')
        patch.set_linewidth('1.5')
    for patch, color in zip(bplot['whiskers'], colors):
        patch.set_color('green') 
        patch.set_linestyle('--')
        patch.set_linewidth('1.5')
    
    # adding horizontal grid lines
    ##for ax in axs:
    #axs.yaxis.grid(True)
    #axs.set_xticks([y + 1 for y in range(len(all_data))])
    #axs.set_xlabel('Four separate samples', fontsize=14)
    #axs.set_ylabel('Observed values')

    # add x-tick labels
    #plt.setp(axs, xticks=[y + 1 for y in range(len(all_data))],
    #        xticklabels=labels)

    #print(mdatadf)
    #top = 1.2
    #bottom = -0.01
    #axs.set_ylim(bottom, top)
    axs.set_xticklabels(labels,
                    rotation=0, fontsize=12) #rotation=45
    
    mdatadf = pd.DataFrame.from_dict(mdatadf)
    '''
    if orient == 'v':
        labels = mdatadf[xl].unique()
        Nlabel = len(labels)
    else:
        labels = mdatadf[yl].unique()
        Nlabel = len(labels)
    '''

    xl="context"
    yl="meth"
    if plotpoint == 'Y':
        ax = sns.swarmplot(data=mdatadf, x=xl, y=yl,
                   alpha=0.5, size=pointsize, #color='grey', size=4,        #颜色，透明度，大小
                   #linewidth=0.5, #edgecolor='gray',        #边线宽度，边线颜色
                   ax=axs,
                   )
    #plt.ylim([-0.01,1.2])


'''
violin plot
'''
def plotviolinplot(fig, axs, mdata, mdatadf, vert, showmeans, showmedians, 
              fontsize, bw_method, showextrema, bxlw, inner, plotpoint, pointsize):

    #fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))

    # Fixing random state for reproducibility
    #np.random.seed(19680801)

    # generate some random test data
    #all_data = [np.random.normal(0, std, 100) for std in range(6, 10)]
    all_data = list(mdata.values())
    labels = list(mdata.keys())
    #print(all_data)

    # plot violin plot
    pos = np.arange(1, len(labels)+1)
    mdatadf = pd.DataFrame.from_dict(mdatadf)
    print(mdatadf, "xxxx")
    xl="context"
    yl="meth"
    #sns.violinplot(data=mdatadf, x=xl, y=yl, split=True, cut=0)
    g=sns.catplot(data=mdatadf, kind='boxen', x=xl, y=yl, 
                      )
    #plt.title("DNA methylation level")
    plt.xlabel("Samples")
    plt.ylabel("DNA methylation level")
    return
    #vplot = axs.violinplot(all_data,
    #                showmeans=showmeans,
    #                showmedians=showmedians,
    #                vert=vert,  # vertical box alignment
    #                showextrema=showextrema, # 
    #                bw_method=bw_method, #str, scalar or callable, optional, This can be 'scott', 'silverman', scott
    #                widths=0.5, #array-like, default: 0.5, Either a scalar or a vector that sets the maximal width of each violin
    #                positions=pos,
    #                )

    #if inner == 'box':
    #    for i,d in enumerate(all_data):
    #        #print(i,d)
    #        min_value,quantile1, median, quantile3, max_value = np.percentile(d, [0, 25, 50, 75, 100])
    #        length = (max_value - min_value)*0.005
    #        print(median)

    #        axs.vlines(i+1, median, median+length, lw=bxlw, zorder=4, color="red", linestyle='-')
    #        ##axs[0].scatter(i+1, median, color='white',zorder=4)
    #        #axs.vlines(i+1,quantile1, quantile3, lw=bxlw, zorder=3, linestyle='-')
    #        #axs.vlines(i+1,min_value, max_value, zorder=2, linestyle='-')
    #        #axs.vlines(i+1, max_value, max_value+length, lw=bxlw, zorder=5, color="red", linestyle='-')
        
    axs.set_title('DNA methylation violin plot')
    
    #colors = ['red', 'lightblue', 'lightgreen', 'red']
    ## add fill color
    #for patch, color in zip(vplot['bodies'], colors):
    #    patch.set_facecolor(color)
    #    patch.set_edgecolor('red') ## 边框
    #    patch.set_linestyle('--') #'-' or 'solid', '--' or 'dashed', '-.' or 'dashdot', ':' or 'dotted', ''
    #    patch.set_linewidth(1.5)

    # adding horizontal grid lines
    #for ax in axs:
    axs.yaxis.grid(True)
    axs.set_xticks([y + 1 for y in range(len(all_data))])
    axs.set_xlabel('Four separate samples')
    axs.set_ylabel('Observed values')

    ## add x-tick labels
    #plt.setp(axs, xticks=[y + 1 for y in range(len(all_data))],
    #        xticklabels=labels)
    ##plt.show()

    if plotpoint == 'Y':
        mdatadf = pd.DataFrame.from_dict(mdatadf)
        xl="context"
        yl="meth"
        ax = sns.swarmplot(data=mdatadf, x=xl, y=yl,
                   alpha=0.5, size=pointsize, #color='grey', alpha=0.5, size=4,        #颜色，透明度，大小
                   #linewidth=0.5, #edgecolor='gray',        #边线宽度，边线颜色
                   ax=axs,
                   )
    


from random import sample

if __name__ == '__main__':
    #plotboxplot()
    args = getArgs().parse_args()
    mydpi = args.dpi
    image_format=args.image_format
    filesplit = args.outfilename.split('.')
    outfileprefix = filesplit[0]
    for i in range(1, len(filesplit)-1):
        outfileprefix = outfileprefix + "." + filesplit[i]
    if image_format == '':
        image_format = filesplit[-1]
    legend=args.legend #是否plot标注
    figextend = 0
    if legend == 11:
        figextend = 1
    #lastlegend = args.lastlegend
    label=args.label
    
    if args.mrfile=='' and args.coverfile == '' and args.statsfile == '':
        print("should predifined input mr files!")
        sys.exit(0)

    if args.coverfile != '':
        coverfs = args.coverfile
        Nfiles = len(coverfs)
        NCcovergae=[]
        xdata=[]
        for cfile in coverfs:
            readcoverfile(cfile, xdata, NCcovergae)
        
        print(xdata, NCcovergae, len(NCcovergae))
        outfilename = outfileprefix + ".cover." + image_format
        plotline(xdata, NCcovergae, outfilename, mydpi, image_format, legend, label)

    if args.statsfile != '':
        statsf = args.statsfile
        Nfiles = len(statsf)
        Sydaya=[]
        Slydata=[]
        Sxdata=[]
        for cfile in statsf:
            readstatsfile(cfile, Sxdata, Sydaya, Slydata)

        if len(label) < Nfiles:
            print("Label number < files, so use filename as labels")
            label = []
            for cfile in statsf:
                label.append(cfile)
        elif len(label) > Nfiles:
            label = label[:Nfiles]

        outfilename = outfileprefix + ".stats." + image_format
        #category_names = [[x] for x in Sxdata]
        category_names = ['[0-0.2)', '[0.2-0.4)',
                              '[0.4-0.6)', '[0.6-0.8)', '[0.8-1]']
        #results = {'m1': Sydaya[0], 'm2': Sydaya[1]}
        results={}
        for i in range(0, len(Sydaya)):
            results[label[i]] = Sydaya[i]
        print('results', results, category_names, Sxdata, Slydata)
        plotline11(Sxdata, Slydata, outfilename, mydpi, image_format, legend, label)
        if len(Slydata[0]) == 5:
            outfilename = outfileprefix + ".stats.bar." + image_format
            outfilename2 = outfileprefix + ".stats.barh." + image_format
            plotbarh(results, category_names, outfilename2, mydpi, image_format, legend, legendsize)
            # 计算每个数组的总和
            sums = [sum(x) for x in Sydaya]
            
            # 对每个数组进行除法运算
            for i in range(len(Sydaya)):
                for j in range(len(Sydaya[i])):
                    Sydaya[i][j] /= sums[i]

            plotbar(label, Sydaya, Slydata, label, outfilename, mydpi, image_format, legend, legendsize)

    if args.mrfile == '':
        exit()
    methlevel = {}
    methleveldf = {}
    mefiles = args.mrfile
    Nmrfile = len(mefiles)
    if len(label) < Nmrfile:
        print("Label number < files, so use filename as labels")
        label = []
        for cfile in mefiles:
            label.append(cfile)
    elif len(label) > Nmrfile:
        label = label[:Nmrfile]

    chr_select = args.chr_select
    i=0
    for mefile in mefiles:
        newframe = readmethfile(methlevel, methleveldf, mefile, label[i], chr_select, len(label), i)
        i+=1
    ##methleveldf = frame
    print(newframe[0].values)
    newmethlevel = {}
    for i in range(len(label)):
        for x in newframe[i].values:
            print(x)
            newmethlevel.setdefault(label[i],[]).append(float(x))
    #print(chr_select, methlevel)
    colors = sample(all_colors, Nmrfile)
    #print(colors)
    notch = False
    vert = True
    showbox = True
    showfliers = False
    showmeans = True
    showmedians = True
    boxprops = {'color':'orangered','facecolor':'pink'}
    medianprops = {}
    meanprops = {}
    capprops = {}
    whiskerprops = {}
    fontsize = 20
    # for violin
    bw_method = 'scott'
    showextrema = True
    bxlw = 9
    inner = 'box---'
    #plotboxplotsns(methlevel)
    #plotboxsns(methlevel, "meth", "context", 'h')
    #plotboxsns(methlevel, "context", "meth", 'v')
    #plotviolionsns(methlevel, "meth", "context", 'h')
    #plotviolionsns(methlevel, "context", "meth", 'v')
    plotmode = args.plotmode
    plotpoint = args.plotpoint
    pointsize = args.pointsize

    outfilename = outfileprefix + ".boxp." + image_format
    if plotmode == "all":
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
        plotboxplot(fig, axs[0], methlevel, methleveldf, notch, vert, showbox, showfliers, showmeans, showmedians, boxprops, medianprops,
                      meanprops, capprops, whiskerprops, fontsize, plotpoint, pointsize)
        plotviolinplot(fig, axs[1], methlevel, methleveldf, vert, showmeans, showmedians,
                      fontsize, bw_method, showextrema, bxlw, inner, plotpoint, pointsize)
    else:
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
        if plotmode == "boxplot":
            plotboxplot(fig, axs, methlevel, methleveldf, notch, vert, showbox, showfliers, showmeans, showmedians, boxprops, medianprops,
                           meanprops, capprops, whiskerprops, fontsize, plotpoint, pointsize)
        elif plotmode == "boxnplot":
            plotviolinplot(fig, axs, methlevel, methleveldf, vert, showmeans, showmedians,
                          fontsize, bw_method, showextrema, bxlw, inner, plotpoint, pointsize)
    fig.tight_layout()
    plt.savefig(outfilename, dpi=mydpi, format=image_format)
    #plt.show()
    
    if len(methlevel) > 1:
        outfilename = outfileprefix + ".corplot1." + image_format
        corplot(newmethlevel, outfilename, mydpi, image_format)
        #pointcorplot(methlevel)
        outfilename = outfileprefix + ".corplot2." + image_format
        pointcorplot(newmethlevel, outfilename, mydpi, image_format)
    
