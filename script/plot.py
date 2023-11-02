import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import os
import sys
import argparse

from util import _load_snp_matrix

def _load_data(fname, header = 1):
    # load matrix file 
    try:
        matrix = np.genfromtxt(fname, skip_header = header, dtype = 'float')
    except IOError:
        print('cannot open {} for reading'.format(fname), file = sys.stderr)
        raise SystemExit
    if len(matrix[:, 1]) != len(matrix[1, :]):
        print('{} {}'.format(len(matrix[:, 1]), len(matrix[1, :])), file = sys.stderr)
        print('matrix file {} is incorrect'.format(fname), file = sys.stderr)
        raise SystemExit
    row, col = np.diag_indices_from(matrix)
    for i in range(1, len(row), 2):
        matrix[i - 1, i - 1] = 0
        matrix[i - 1, i] = 0
        matrix[i, i - 1] = 0
        matrix[i, i] = 0
    for i in range(1, len(row), 2):
        for j in range(1, len(col), 2):
            total = matrix[i - 1, j - 1] + matrix[i - 1, j] + matrix[i, j - 1] + matrix[i, j] + 0.000001
            matrix[i - 1, j - 1] /= total
            matrix[i - 1, j] /= total
            matrix[i, j - 1] /= total
            matrix[i, j] /= total
    return matrix

def heat_map_plot(matrix, name, ofname):
    # matrix = _load_data(fname, header)
    length = len(matrix)
    plt.figure(figsize = (6, 6), facecolor = 'w', edgecolor = 'w')
    with np.errstate(divide = 'ignore'):
        ax = plt.imshow(matrix, cmap = plt.get_cmap('YlOrRd'), interpolation = 'nearest', aspect = 'equal', extent = (0, length, 0, length), origin = 'lower')
    ticks = [i for i in range(0, len(matrix), 2)]
    plt.xticks(ticks, [])
    plt.yticks(ticks, [])
    plt.grid(color = 'black', linewidth = 0.5, alpha = 0.5)
    plt.tick_params(direction = 'in')
    plt.title(name)
    divider = make_axes_locatable(ax.axes)
    cax = divider.append_axes('bottom', size = '2.5%', pad = 0.3)
    plt.colorbar(ax, cax = cax, ticks = MultipleLocator(0.2), format = '%.1f',  orientation = 'horizontal', extendfrac = 'auto', spacing = 'uniform')
    plt.savefig(ofname, dpi = 600)
    plt.clf()
    plt.close()

def heat_map_all_plot(mfname, dir):
    os.makedirs(dir, exist_ok = True)
    for header, matrix in _load_snp_matrix(mfname):
        ctg = header[0].split(':')[0].split('_')[0]
        matrix_np = np.array(matrix, dtype = 'float')
        row, col = np.diag_indices_from(matrix_np)
        for i in range(1, len(row), 2):
            matrix_np[i - 1, i - 1] = 0
            matrix_np[i - 1, i] = 0
            matrix_np[i, i - 1] = 0
            matrix_np[i, i] = 0
        for i in range(1, len(row), 2):
            for j in range(1, len(col), 2):
                total = matrix_np[i - 1, j - 1] + matrix_np[i - 1, j] + matrix_np[i, j - 1] + matrix_np[i, j] + 0.000001
                matrix_np[i - 1, j - 1] /= total
                matrix_np[i - 1, j] /= total
                matrix_np[i, j - 1] /= total
                matrix_np[i, j] /= total
        ofname = os.path.join(dir, '{}.png'.format(ctg))
        heat_map_plot(matrix_np, ctg, ofname)

def cov_plot(cfname, windows, dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    files = []
    for fpath, dirs, fs in os.walk(dir):
        for f in fs:
            if cfname not in f: continue
            filename = os.path.join(fpath, f)
            if filename.endswith('.cov'):
                files.append(filename)

    for file in files:
        # plot coverage figure for every contig
        # windows = 1
        cov = np.loadtxt(file, dtype = np.int)
        if not np.any(cov): continue
        fname = os.path.basename(file).replace('.cov', '.png')
        fname = os.path.join(dir, fname)
        cov = np.convolve(cov, np.ones(windows), 'valid') / windows
        plt.plot(cov)
        plt.savefig(fname)
        plt.clf()

def cov_plot2(cfname, windows, output):
    cov_diphase = np.loadtxt(cfname[0], dtype = np.int)
    if not np.any(cov_diphase): return
    cov_diphase = np.convolve(cov_diphase, np.ones(windows), 'valid') / windows
    plt.plot(cov_diphase, label = 'Diphase')
    cov_falcon = np.loadtxt(cfname[1], dtype = np.int)
    if not np.any(cov_falcon): return
    cov_falcon = np.convolve(cov_falcon, np.ones(windows), 'valid') / windows
    plt.plot(cov_falcon, label = 'FALCON-Phase')
    plt.legend()
    plt.savefig(output)

def heatmap(args):
    matrix = _load_data(args.fname, args.header)
    heat_map_plot(matrix, args.name, args.out, args.header)

def heatmap_all(args):
    heat_map_all_plot(args.matrix, args.dir)

def cov(args):
    cov_plot2(args.cov, args.windows, args.dir)

def main():
    parser = argparse.ArgumentParser('plot')
    subparser = parser.add_subparsers(title = 'subcommand')

    parser_heatmap = subparser.add_parser('heatmap', help = 'plot heatmap')
    parser_heatmap.add_argument('-f', '--fname', dest = 'fname', default = None, metavar = 'path', help = 'matrix file')
    parser_heatmap.add_argument('-n', '--name', dest = 'name', default = None, metavar = 'str', help = 'title of figure')
    parser_heatmap.add_argument('-o', '--out', dest = 'out', default = 'phasing', metavar = 'str', help = 'output file name')
    parser_heatmap.add_argument('-hd', '--header', dest = 'header', default = 1, metavar = 'int', help = 'number of line(s) should be ignore in the header of the matrix file (default 1)')
    parser_heatmap.set_defaults(func = heatmap)

    parser_heatmap_all = subparser.add_parser('heatmap_all', help = 'plot heatmap for all contigs')
    parser_heatmap_all.add_argument('-m', '--matrix', dest = 'matrix', default = None, metavar = 'path', help = 'matrix file')
    parser_heatmap_all.add_argument('-d', '--dir', dest = 'dir', default = None, metavar = 'path', help = 'directory to save figures')
    parser_heatmap_all.set_defaults(func = heatmap_all)

    parser_cov = subparser.add_parser('cov', help = 'coverage plot of contigs based on Hi-C mapping')
    parser_cov.add_argument('-c', '--cov', dest = 'cov', default = None, nargs = '+', metavar = 'path', help = 'coverage information file name')
    parser_cov.add_argument('-w', '--windows', dest = 'windows', type = int, default = 1, metavar = 'int', help = 'windows size for smoothing coverage (default 1)')
    parser_cov.add_argument('-d', '--dir', dest = 'dir', default = None, metavar = 'path', help = 'directory to save figures')
    parser_cov.set_defaults(func = cov)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()