#!/usr/bin/env python
# -*- coding: utf-8 -*-

## anno. this script modified based on plotCorrelation in deeptools
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import os
from correlation import Correlation

old_settings = np.seterr(all='ignore')

def writableFile(string):
    """
    Simple function that tests if a given path is writable
    """
    try:
        open(string, 'w').close()
        os.remove(string)
    except:
        msg = "{} file can't be opened for writing".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string

def parse_arguments(args=None):
    basic_args = plot_correlation_args()
    heatmap_parser = heatmap_options()
    scatter_parser = scatterplot_options()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Tool for the analysis and visualization of sample correlations based on the output of DMtools. 
Pearson or Spearman methods are available to compute correlation coefficients. 
Results can be saved as multiple scatter plots depicting the pairwise correlations or as a clustered heatmap,
where the colors represent the correlation coefficients and the clusters are constructed using complete linkage.
Optionally, the values can be saved as tables, too.

detailed help:

  DMplotcor -h

""",
        epilog='example usages:\n'
               'DMplotcor -i results_file -p heatmap -c pearson -o heatmap.pdf\n'
               'or\n'
               'DMplotcor -f sample1.bins.file sample2.bins.file sample3.bins.file -p heatmap -c pearson -o heatmap.pdf\n'
               ' \n\n',
        parents=[basic_args, heatmap_parser, scatter_parser])

    return parser


def plot_correlation_args():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    parser.add_argument('-f', '--binFiles', nargs='+',
                            help='the list of input DNA methylation files to merge and plot, must provide -f or -i')

    required.add_argument('--corMatrix', '-i',
                          metavar='FILE',
                          help='Matrix of DNA methylation values, must provide -f or -i')

    required.add_argument('--corMethod', '-m',
                          help="Correlation method.",
                          choices=['spearman', 'pearson'],
                          default='spearman',
                          required=True)

    required.add_argument('--whatToPlot', '-p',
                          help="Choose between a heatmap or pairwise scatter plots",
                          choices=['heatmap', 'scatterplot'],
                          default='heatmap',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--plotFile', '-o',
                          help='File to save the heatmap to. The file extension determines the format, '
                          'so heatmap.pdf will save the heatmap in PDF format. '
                          'The available formats are: .png, '
                          '.eps, .pdf and .svg.',
                          type=writableFile,
                          metavar='FILE')

    required.add_argument('--context', '-c',
                          help='DNA methylation context for calculate, All equal to CG+CHG+CHH, '
                          'C means average of CG/CHG/CHH',
                          choices=['All', 'C', 'CG', 'CHG', 'CHH'],
                          default='heatmap',
                          required=True)

    optional.add_argument('--skipZeros',
                          help='By setting this option, genomic regions '
                          'that have zero or missing (nan) values in all samples '
                          'are excluded.',
                          action='store_true',
                          required=False)

    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                          'file names. '
                          'Multiple labels have to be separated by spaces, e.g. '
                          '--labels sample1 sample2 sample3',
                          nargs='+')

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title. (Default: %(default)s)',
                          default='')

    optional.add_argument('--plotFileFormat', '-pf',
                          metavar='FILETYPE',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. The available options are: png, '
                          'eps, pdf and svg.',
                          choices=['png', 'pdf', 'svg', 'eps', 'plotly'])

#    optional.add_argument(
#        '--removeOutliers',
#        help='If set, bins with very large counts are removed. ',
#        action='store_true')

#    optional.add_argument('--version', action='version',
#                          version='%(prog)s {}'.format(__version__))

    group = parser.add_argument_group('Output optional options')

    group.add_argument('--outFileCorMatrix',
                       help='Save matrix with pairwise correlation values to a tab-separated file.',
                       metavar='FILE',
                       type=writableFile)

    return parser


def scatterplot_options():
    """
    Options specific for creating the scatter plot
    """
    parser = argparse.ArgumentParser(add_help=False)
    scatter_opts = parser.add_argument_group('Scatter plot options')

    scatter_opts.add_argument('--xRange',
                              help='The X axis range. The default scales these such that the full range of dots is displayed.',
                              type=float,
                              nargs=2,
                              default=None)

    scatter_opts.add_argument('--yRange',
                              help='The Y axis range. The default scales these such that the full range of dots is displayed.',
                              type=float,
                              nargs=2,
                              default=None)

#    scatter_opts.add_argument('--log1p',
#                              help='Plot the natural log of the scatter plot after adding 1. Note that this is ONLY for plotting, the correlation is unaffected.',
#                              action='store_true')

    return parser


def heatmap_options():
    """
    Options for generating the correlation heatmap
    """
    parser = argparse.ArgumentParser(add_help=False)
    heatmap = parser.add_argument_group('Heatmap options')

    heatmap.add_argument('--plotHeight',
                         help='Plot height in cm. (Default: %(default)s)',
                         type=float,
                         default=9.5)

    heatmap.add_argument('--plotWidth',
                         help='Plot width in cm. The minimum value is 1 cm. (Default: %(default)s)',
                         type=float,
                         default=11)

    heatmap.add_argument('--zMin', '-min',
                         default=None,
                         help='Minimum value for the heatmap intensities. '
                              'If not specified, the value is set automatically',
                         type=float)

    heatmap.add_argument('--zMax', '-max',
                         default=None,
                         help='Maximum value for the heatmap intensities.'
                              'If not specified, the value is set automatically',
                         type=float)

    heatmap.add_argument(
        '--colorMap', default='jet',
        metavar='',
        help='Color map to use for the heatmap. Available values can be '
             'seen here: '
             'http://matplotlib.org/examples/color/colormaps_reference.html')

    heatmap.add_argument('--plotNumbers',
                         help='If set, then the correlation number is plotted '
                         'on top of the heatmap. This option is only valid when plotting a heatmap.',
                         action='store_true',
                         required=False)

    return parser


#def main(args=None):
if __name__ == "__main__":
    args = parse_arguments().parse_args()

    if args.plotFile is None and args.outFileCorMatrix is None:
        sys.exit("At least one of --plotFile and --outFileCorMatrix must be specified!\n")

    if args.corMatrix is None and args.binFiles is None:
        sys.exit("At least one of -i/--corMatrix and -f/--binFiles must be specified!\n")

    args.removeOutliers = False
    args.log1p = False
    corr = Correlation(args.corMatrix, args.binFiles,
                       args.corMethod,
                       labels=args.labels,
                       remove_outliers=args.removeOutliers,
                       skip_zeros=args.skipZeros,
                       context=args.context)

    if args.corMethod == 'pearson':
        # test if there are outliers and write a message recommending the removal
        if len(corr.get_outlier_indices(np.asarray(corr.matrix).flatten())) > 0:
            if args.removeOutliers:
                sys.stderr.write("\nOutliers were detected in the data. They "
                                 "will be removed to avoid bias "
                                 "in the pearson correlation.\n")

            else:
                sys.stderr.write("\nOutliers were detected in the data. Consider "
                                 "using the --removeOutliers parameter to avoid a bias "
                                 "in the pearson correlation.\n")

    if args.colorMap:
        try:
            plt.get_cmap(args.colorMap)
        except ValueError as error:
            sys.stderr.write(
                "A problem was found. Message: {}\n".format(error))
            exit()

    if args.plotFile is not None:
        if args.whatToPlot == 'scatterplot':
            corr.plot_scatter(args.plotFile,
                              plot_title=args.plotTitle,
                              image_format=args.plotFileFormat,
                              xRange=args.xRange,
                              yRange=args.yRange,
                              log1p=args.log1p)
        else:
            corr.plot_correlation(args.plotFile,
                                  vmax=args.zMax,
                                  vmin=args.zMin,
                                  colormap=args.colorMap,
                                  plot_title=args.plotTitle,
                                  image_format=args.plotFileFormat,
                                  plot_numbers=args.plotNumbers,
                                  plotWidth=args.plotWidth,
                                  plotHeight=args.plotHeight)

    if args.outFileCorMatrix:
        o = open(args.outFileCorMatrix, "w")
        o.write("#DMplotcor --outFileCorMatrix\n")
        corr.save_corr_matrix(o)
        o.close()

