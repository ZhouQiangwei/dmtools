# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 10:29:50 2018

@author: qwzhou
=======================================
plot line and dash
=======================================
This is a script for profile across the TSS/TES
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.backends.backend_pdf import PdfPages
import argparse

plt.switch_backend('agg')

y = []
z = []
k = []
averm = []

def readchromfile(filename, Nchr):
    chr=''
    with open(filename, 'r') as fig:
        for line in fig:
            data = line.split()
            ## chrom
            #Chr1    0       +       0.022330        C
            ## chrom.cover.c
            #Chr1    0       +       15685   702419  Chr1:0-100000
            chrom, pos, strand, meth, context = data
            if chrom not in Nchr:
                Nchr[chrom] = 0
            Nchr[chrom] = Nchr[chrom]+1
            

######################################################
def find_martrix_max_value(data_matrix):  
    new_data=[]  
    for i in range(len(data_matrix)):  
        new_data.append(max(data_matrix[i]))
    return max(new_data)

skipyaxis = False
plotonly = False
Nfigure=3
legendsize=10
#######################################################
def plotline( plotx, ploty, title, label, nsample, legend, filename, scale_ls, index_ls, ylabel, plotn, mycolor,ymin,ymax):
    prosamp = 0
    if plotonly:
        ax = fig.subplots()
    else:
        ax = fig.add_subplot(1,Nfigure,plotn) #add_subplot(nrows, ncols, index)
    while prosamp < nsample:
        ploty[prosamp] = [i*percentage for i in ploty[prosamp]]
        #for i,item in enumerate(y[prosamp]):
        #    if item >6:
        #        y[prosamp][i]=6
        if(len(mycolor)>=nsample):
            ax.plot(plotx, ploty[prosamp], label=label[prosamp], color=mycolor[prosamp])
        else:
            #print(prosamp, label,  label[prosamp])
            ax.plot(plotx, ploty[prosamp], label=label[prosamp]) 
        prosamp = prosamp +1
    #dashes = [10, 5, 100, 5]
    #line1.set_dashes(dashes)   #   dash line
    # Remove the plot frame lines. They are unnecessary here.
    #ax.spines['top'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter('{:.0f}'.format))
    #ax.yaxis.set_major_formatter(plt.FuncFormatter('{:.1f}%'.format))
    #plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                labelbottom=True, left=True, right=False, labelleft=True)
    #ax.axes.get_xaxis().set_visible(False)
    if legend < 11:
        plt.legend(loc=legend, prop={'size': legendsize}, frameon=False)
    elif legend == 11:
        plt.legend(bbox_to_anchor=(1.05, 0.05), loc=3, borderaxespad=0, prop={'size': legendsize}, frameon=False)
    
    #plt.axvline(x=xlen/2-1, ls="--", color='black')
    
    #plt.axhline(y=0, xmin=0.05, xmax=0.5, linewidth=8, color='gray')
    #plt.axhline(y=0, xmin=0.5, xmax=0.505, linewidth=8, color='k' )
    #plt.axhline(y=0, xmin=0.505, xmax=0.95, linewidth=8, color='gray') 
    
    #print(scale_ls, index_ls, len(plotx))
    #index_ls = ['-200bp','Start', "+200bp"]
    #plt.xticks(scale_ls2,index_ls,color='k', size=12)
    ax.set_title(title,size=12)
    ax.set_ylabel(ylabel,size=12)
    #ax.set_ylabel('Methylation Level',size=15)
    if skipyaxis:
        miny = 0
        maxy = 1
        maxy = find_martrix_max_value(ploty) * 1.1
    else:
        miny = ymin
        maxy = ymax
    
    ax.set_ylim(miny, maxy)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=45)
    #plt.show()
    
    #filename2=filename + ".png"
    #plt.savefig(filename2, bbox_inches='tight')


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
                        help="DNA AverMethylevel files, seperate by space. eg. wildtype.chrommethy.txt", 
                        nargs='+')
    parser.add_argument("-l", "--label", default='wildtype', 
                        help="Labels of samples, sperate by space. eg. -l widetype", 
                        nargs='+')
    parser.add_argument('--outFileName', '-o',
                       help='Output file name.',
                       metavar='FILENAME',
                       default='myMeth.pdf',
                       type=writableFile,
                       required=True)
    parser.add_argument("-yl", "--ylabel",
                        default="",
                        #nargs='+',
                        help='y-axis label')
    parser.add_argument( '-t', '--title',
                        default=["C meth", "CG meth", "CHG meth", "CHH meth"],
                        help='Title of the plot, to be printed on top of '
                        'the generated image. Leave blank for no title.',
                        nargs='+')
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
    parser.add_argument('--color',
                        default=["red"],
                        nargs='+',
                        help='List of colors to use, should same as the number of samples,'
                        'Color names and html hex strings (e.g., #eeff22) are accepted. The color names should '
                        'be space separated. For example, --color red blue green ')
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
    parser.add_argument('--lastlegend',
                        default=True,
                        help='Only show the last figure\'s legend.'
                        )
    parser.add_argument('--legendsize',
                        default=10,
                        help='the text size of the legend.',
                        type=int
                        )
    parser.add_argument('--context',
                        default="ALL",
                        #choices=["ALL", "C", "CG", "CHG","CHH"],
                        help='List of colors to use, should same as the number of samples,'
                        'Color names and html hex strings (e.g., #eeff22) are accepted. The color names should '
                        'be space separated. For example, --color red blue green ')
    parser.add_argument("-ft", "--image_format",
                        default="pdf",
                        help="The file format, e.g. 'png', 'pdf', 'svg', ... "
                        "The behavior when this is unset is documented under fname.")
    parser.add_argument("--dpi",
                        default=200,
                        help="Set the DPI to save the figure. default: 200",
                        type=int)
    parser.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
                        

    return parser


if __name__ == '__main__':
    #args = vars(ap.parse_args())
    args = getArgs().parse_args()
    #methpoint(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    
    if(args.mrfile==''):
        print("should predifined input mr files!")
        exit()
    
    allposdata = {}
    allmedata = {}
    label=args.label
    ylabel=args.ylabel
    Nmrfile = len(args.mrfile)
    if len(label) < Nmrfile:
        print("Label number < files, so use filename as labels")
        label = []
        for cfile in mefiles:
            label.append(cfile)
    elif len(label) > Nmrfile:
        label = label[:Nmrfile]

    Nchrs = {}
    for x in range(0, len(args.mrfile)):
        readchromfile(args.mrfile[x], Nchrs)

    Nchr=0
    for x in Nchrs:
        if Nchrs[x] > 20:
            Nchr+=1
    print(Nchrs, Nchr)

    filename=args.outFileName
    nsample=len(args.mrfile)
    legend=args.legend #是否plot标注
    figextend = 0
    if legend == 11:
        figextend = 1
    lastlegend = args.lastlegend

    image_format=args.image_format
    legendsize=args.legendsize
    
    title=args.title
    mycolor=args.color
    plotn = 1
    ymin=args.yMin
    ymax=args.yMax
    context=args.context

    newcontexts = context.split(',')
    contextlen = len(newcontexts)

    mydpi=args.dpi

    Nfigure=1
    plotonly = False
    fig = plt.figure(figsize=((4+figextend)*Nchr,3))
    Nfigure=4
    
    if(context=="ALL"):
        if(len(ymax)<4):
            skipyaxis = True
        while(len(ymax)<4):
            ymax.append(1.0)
        while(len(ymin)<4):
            ymin.append(0.0)
    elif contextlen > 1:
        if(len(ymax)<contextlen):
            skipyaxis = True
        while(len(ymax)<contextlen):
            ymax.append(1.0)
        while(len(ymin)<contextlen):
            ymin.append(0.0)
    else:
        if(len(ymax)<1):
            skipyaxis = True
        while(len(ymax)<1):
            ymax.append(1.0)
        while(len(ymin)<1):
            ymin.append(0.0)

    for chr in allmedata:
        print(chr, allmedata[chr], allposdata[chr])
        plotline( plotx, ploty, title, label, nsample, legend, filename, scale_ls, index_ls, ylabel, plotn, mycolor, ymin, ymax)

    plt.subplots_adjust(wspace=0.05, hspace=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=mydpi, format=image_format)
    plt.close()
