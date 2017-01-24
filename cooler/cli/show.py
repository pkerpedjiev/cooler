# -*- coding: utf-8 -*-
from __future__ import division, print_function
import sys

import numpy as np
import h5py

import click
from . import cli
from ..api import Cooler
from .. import util


MAX_MATRIX_SIZE_FILE = int(1e8)
MAX_MATRIX_SIZE_INTERACTIVE = int(1e7)


fallList = ((255, 255, 255), (255, 255, 204),
     (255, 237, 160), (254, 217, 118),
     (254, 178, 76), (253, 141, 60),
     (252, 78, 42), (227, 26, 28),
     (189, 0, 38), (128, 0, 38), (0, 0, 0))

bluesList = ((255, 255, 255), (180, 204, 225),
             (116, 169, 207), (54, 144, 192),
             (5, 112, 176), (4, 87, 135),
             (3, 65, 100), (2, 40, 66),
             (1, 20, 30), (0, 0, 0))


acidBluesList = ((255, 255, 255), (162, 192, 222),
                 (140, 137, 187), (140, 87, 167),
                 (140, 45, 143), (120, 20, 120),
                 (90, 15, 90), (60, 10, 60),
                 (30, 5, 30), (0, 0, 0))

nMethList = ((236, 250, 255), (148, 189, 217),
             (118, 169, 68), (131, 111, 43), (122, 47, 25),
             (41, 0, 20))

def registerList(mylist, name):
    import matplotlib as mpl
    mymap = listToColormap(mylist, name)
    mymapR = listToColormap(mylist[::-1], name + "_r")
    mpl.cm.register_cmap(name, mymap)
    mpl.cm.register_cmap(name + "_r", mymapR)

def listToColormap(colorList, cmapName=None):
    import matplotlib as mpl
    colorList = np.array(colorList)
    if colorList.min() < 0:
        raise ValueError("Colors should be 0 to 1, or 0 to 255")
    if colorList.max() > 1.:
        if colorList.max() > 255:
            raise ValueError("Colors should be 0 to 1 or 0 to 255")
        else:
            colorList = colorList / 255.
    return mpl.colors.LinearSegmentedColormap.from_list(cmapName, colorList, 256)

def gridspec_inches(
    wcols, 
    hrows, 
    fig_kwargs={}):

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    fig_height_inches = (
        sum(hrows) 
        )

    fig_width_inches = (
        sum(wcols)
        )

    fig=plt.figure(
        figsize=(fig_width_inches,fig_height_inches),
        subplotpars=mpl.figure.SubplotParams(
        left=0,
        right=1,
        bottom=0,
        top=1,
        wspace =0,
        hspace = 0.0),
        frameon=False,
        **fig_kwargs)
    fig.set_size_inches(fig_width_inches,fig_height_inches,forward=True)

    gs = mpl.gridspec.GridSpec(
        len(hrows),
        len(wcols),
        left=0,
        right=1,
        top=1,
        bottom=0,
        wspace=0,
        hspace=0,
        width_ratios=wcols,
        height_ratios=hrows
        )

    return fig, gs

def get_matrix_size(c, row_region, col_region):
    nrows = c.extent(row_region)[1] - c.extent(row_region)[0]
    ncols = c.extent(col_region)[1] - c.extent(col_region)[0]
    return ncols * nrows

def load_matrix(c, row_region, col_region, balanced, scale):
    mat = (c.matrix(balance=balanced)
            .fetch(row_region, col_region)
            .toarray())

    if scale == 'log2':
        mat = np.log2(mat)
    elif scale == 'log10':
        mat = np.log10(mat)

    return mat

def load_track(path, region):
    import pyBigWig
    region_parsed = util.parse_region_string(region)
    bw = pyBigWig.open(path)
    track = bw.values(*region_parsed)
    xs = np.arange(region_parsed[1],region_parsed[2])
    return xs, track, region_parsed

def interactive(plotstate, c, row_chrom, col_chrom, balanced, scale):
    import matplotlib.pyplot as plt
    # The code is heavily insired by 
    # https://gist.github.com/mdboom/048aa35df685fe694330764894f0e40a
    
    def get_extent(ax):
        xstart, ystart, xdelta, ydelta = ax.viewLim.bounds
        xend = xstart + xdelta
        yend = ystart + ydelta
        return xstart, xend, ystart, yend

    def round_trim_extent(extent, binsize, row_chrom_len, col_chrom_len):
        xstart = int(np.floor(extent[0] / binsize) * binsize)
        xend = int(np.ceil(extent[1] / binsize) * binsize)
        ystart = int(np.floor(extent[3] / binsize) * binsize)
        yend = int(np.ceil(extent[2] / binsize) * binsize)
        xstart = max(0, xstart)
        ystart = max(0, ystart)
        # For now, don't let users to request the last bin, b/c its end
        # lies outside of the genome
        xend = min(xend, int(np.floor(col_chrom_len / binsize) * binsize))
        yend = min(yend, int(np.floor(row_chrom_len / binsize) * binsize))
        return xstart, xend, yend, ystart

    def move_data(event):
        for ax in plotstate['data_axes']['hms']:
            ax.set_autoscale_on(False)  # Otherwise, infinite loop

        extent = get_extent(event.inaxes)
        extent = round_trim_extent(extent, binsize, row_chrom_len, col_chrom_len)
        if event.inaxes in plotstate['data_axes']['tracks']:
            extent = (
                extent[0], 
                extent[1], 
                plotstate['prev_extent'][2],
                plotstate['prev_extent'][3]
            )

        print(extent)
        if (extent == plotstate['prev_extent']):
            return

        plotstate['prev_extent'] = extent
        new_col_region = col_chrom, int(extent[0]), int(extent[1])
        new_row_region = row_chrom, int(extent[3]), int(extent[2])

        for ax in plotstate['data_axes']['hms']:
            im = ax.images[-1]
            nelem = get_matrix_size(c, new_row_region, new_col_region)
            if nelem  >= MAX_MATRIX_SIZE_INTERACTIVE:
                # requested area too large
                im.set_data(np.ones(1)[:, None] * np.nan)
                
                if not plotstate['placeholders']:
                    box, = plt.plot(
                        [0, col_chrom_len, col_chrom_len, 0, 0, col_chrom_len],
                        [0, row_chrom_len, 0, 0, row_chrom_len, row_chrom_len],
                        c='k',
                        lw=0.5
                    )
                    txt = plt.text(
                        0.5, 0.5,
                        'The requested region is too large\n'
                        'to display at this resolution.',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes
                    )
                    plotstate['placeholders'] = [box, txt]
            else:
                # remove placeholders if any and update
                while plotstate['placeholders']:
                    plotstate['placeholders'].pop().remove()
                
                im.set_data(
                    load_matrix(c, new_row_region, new_col_region, balanced, scale))

            im.set_extent(extent)
            ax.set_xlim(*extent[:2])
            ax.set_ylim(*extent[-2:])
            ax.figure.canvas.draw_idle()

        for track_line in plotstate['track_lines']:
            xs, track, region_parsed = load_track(track_line['path'],
                '{}:{}-{}'.format(*new_col_region))
            track_line['line'].set_data(xs, track)
            track_line['ax'].set_xlim(*extent[:2])


    binsize = c.info['bin-size']
    chromsizes = c.chroms()[:].set_index('name')['length']
    row_chrom_len = chromsizes[row_chrom]
    col_chrom_len = chromsizes[col_chrom]
    plotstate['placeholders'] = []
    plotstate['prev_extent'] = get_extent(plotstate['data_axes']['hms'][0])
    plt.gcf().canvas.mpl_connect('button_release_event', move_data)
    plt.show()


@cli.command()
@click.argument(
    "cooler_file",
    metavar="COOLER_PATH")
@click.argument(
    "range",
    type=str)
@click.option(
    "--range2", "-r2",
    type=str,
    help="The coordinates of a genomic region shown along the column dimension. "
         "If omitted, the column range is the same as the row range. "
         "Use to display asymmetric matrices or trans interactions.")
@click.option(
    "--balanced", "-b",
    is_flag=True,
    default=False,
    help="Show the balanced contact matrix. "
         "If not provided, display the unbalanced counts.")
@click.option(
    "--out", "-o",
    help="Save the image of the contact matrix to a file. "
         "If not specified, the matrix is displayed in an interactive window. "
         "The figure format is deduced from the extension of the file, "
         "the supported formats are png, jpg, svg, pdf, ps and eps.")
@click.option(
    "--dpi",
    type=int,
    help="The DPI of the figure, if saving to a file")
@click.option('--scale', '-s',
    type=click.Choice(['linear', 'log2', 'log10']),
    help="Scale transformation of the colormap: linear, log2 or log10. "
         "Default is log10.",
    default='log10')
@click.option(
    "--force", "-f",
    is_flag=True,
    default=False,
    help="Force display very large matrices (>=10^8 pixels). "
         "Use at your own risk as it may cause performance issues.")
@click.option(
    "--zmin",
    type=float,
    help="The minimal value of the color scale. Units must match those of the colormap scale. "
         "To provide a negative value use a equal sign and quotes, e.g. -zmin='-0.5'")
@click.option(
    "--zmax",
    type=float,
    help="The maximal value of the color scale. Units must match those of the colormap scale. "
         "To provide a negative value use a equal sign and quotes, e.g. -zmax='-0.5'")
@click.option(
    "--cmap",
    #default="YlOrRd",
    default="fall",
    help="The colormap used to display the contact matrix. "
         "See the full list at http://matplotlib.org/examples/color/colormaps_reference.html")

@click.option(
    "--track",
    type=(str, int, str),
    help="Extra track to display",
    multiple=True
    )

def show(cooler_file, 
    range, range2, balanced, out, dpi, scale, 
    force, zmin, zmax, cmap,
    track):
    """
    Display a contact matrix.
    Display a region of a contact matrix stored in a COOL file.

    COOLER_PATH : Path to a COOL file

    RANGE : The coordinates of the genomic region to display, in UCSC notation.
    Example: chr1:10,000,000-11,000,000

    """
    try:
        import matplotlib as mpl
        if out is not None:
            mpl.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("Install matplotlib to use cooler show", file=sys.stderr)
        sys.exit(1)

    registerList(fallList, "fall")
    registerList(bluesList, "blues")
    registerList(acidBluesList, "acidblues")
    registerList(nMethList, "nmeth")

    c = Cooler(cooler_file)
    chromsizes = c.chroms()[:].set_index('name')['length']
    row_region = range
    col_region = row_region if range2 is None else range2
    row_chrom, row_lo, row_hi = util.parse_region(row_region, chromsizes)
    col_chrom, col_lo, col_hi = util.parse_region(col_region, chromsizes)

    if ((get_matrix_size(c, row_region, col_region) >= MAX_MATRIX_SIZE_FILE) 
        and not force):
        print(
            "The matrix of the selected region is too large. "
            "Try using lower resolution, selecting a smaller region, or use "
            "the '--force' flag to override this safety limit.",
            file=sys.stderr)
        sys.exit(1)

    n_track_windows = max([1+t[1] for t in list(track)]) if len(list(track))>0 else 0
    fig, gs = gridspec_inches(
        [1,8,1,1], 
        [0.5,8,0.5] + [1, 0.5]*n_track_windows
    )
    plotstate = {}
    plotstate['data_axes'] = {'hms':[], 'tracks':[]}
    plotstate['track_lines'] = []
    fig.canvas.set_window_title('Contact matrix'.format())
    plt.subplot(gs[1,1])
    plotstate['data_axes']['hms'].append(plt.gca())
    plt.title('')
    hm = plt.imshow(
        load_matrix(c, row_region, col_region, balanced, scale),
        interpolation='none',
        extent=[col_lo, col_hi, row_hi, row_lo],
        vmin=zmin,
        vmax=zmax,
        cmap=cmap)

    # If plotting into a file, plot and quit
    plt.ylabel('{} coordinate'.format(row_chrom))
    plt.xlabel('{} coordinate'.format(col_chrom))

    plt.subplot(gs[1,2])
    plt.axis('off')
    cb = plt.colorbar(hm, aspect=20, fraction=0.8)
    cb.set_label(
        {'linear': 'relative contact frequency',
         'log2'  : 'log 2 ( relative contact frequency )',
         'log10' : 'log 10 ( relative contact frequency )'}[scale])

    if track:
        n_tracks = {}
        for path, i, track_title in track:
            plt.subplot(gs[3+i*2, 1])
            if plt.gca() not in plotstate['data_axes']['tracks']:
                plotstate['data_axes']['tracks'].append(plt.gca())
            n_tracks[i] = n_tracks.get(i, 0) + 1 

            xs, track, region_parsed = load_track(path, range)
            line, = plt.plot(
                xs,
                track,
                label=track_title
            )

            plotstate['track_lines'].append({
                'line':line, 
                'ax':plt.gca(), 
                'path':path, 
                'idx':i, 
                'title':track_title})
            
            plt.legend(loc='upper center', ncol=n_tracks[i])
        
    if out:
        plt.savefig(out, dpi=dpi)
    else:
        interactive(plotstate, c, row_chrom, col_chrom, balanced, scale)
