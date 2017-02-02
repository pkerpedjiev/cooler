# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import six
import os

from scipy.sparse import coo_matrix
import numpy as np
import pandas
import h5py

from .core import RangeSelector1D, RangeSelector2D, TriuReader, query_rect
from .util import parse_region, open_hdf5


def get(h5, lo=0, hi=None, fields=None, convert_enum=True, **kwargs):
    """
    Query a range of rows from a table as a dataframe.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    h5 : ``h5py.Group`` or any dict-like of array-likes
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : str or sequence of str, optional
        Column or list of columns to query. Defaults to all available columns.
        A single string returns a Series instead of a DataFrame.
    convert_enum : bool, optional
        Whether to convert HDF5 enum datasets into ``pandas.Categorical``
        columns instead of plain integer columns. Default is True.
    kwargs : optional
        Options to pass to ``pandas.DataFrame`` or ``pandas.Series``.

    Returns
    -------
    DataFrame or Series

    Notes
    -----
    HDF5 ASCII datasets are converted to Unicode.

    """
    grp = h5
    series = False
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, six.string_types):
        fields = [fields]
        series = True

    data = {}
    for field in fields:
        dset = grp[field]

        if convert_enum:
            dt = h5py.check_dtype(enum=dset.dtype)
        else:
            dt = None

        if dt is not None:
            data[field] = pandas.Categorical.from_codes(
                dset[lo:hi],
                sorted(dt, key=dt.__getitem__),
                ordered=True)
        elif dset.dtype.type == np.string_:
            data[field] = dset[lo:hi].astype('U')
        else:
            data[field] = dset[lo:hi]

    if data and lo is not None:
        index = np.arange(lo, lo + len(next(iter(data.values()))))
    else:
        index = None

    if series:
        return pandas.Series(
            data[fields[0]],
            index=index,
            name=field,
            **kwargs)
    else:
        return pandas.DataFrame(
            data,
            columns=fields,
            index=index,
            **kwargs)


def _region_to_extent(h5, chrom_ids, region, binsize):
    chrom, start, end = region
    cid = chrom_ids[chrom]
    if binsize is not None:
        chrom_offset = h5['indexes']['chrom_offset'][cid]
        yield chrom_offset + int(np.floor(start/binsize))
        yield chrom_offset + int(np.ceil(end/binsize))
    else:
        chrom_lo = h5['indexes']['chrom_offset'][cid]
        chrom_hi = h5['indexes']['chrom_offset'][cid + 1]
        chrom_bins = h5['bins']['start'][chrom_lo:chrom_hi]
        yield chrom_lo + np.searchsorted(chrom_bins, start, 'right') - 1
        yield chrom_lo + np.searchsorted(chrom_bins, end, 'left')


def region_to_offset(h5, chrom_ids, region, binsize=None):
    return next(_region_to_extent(h5, chrom_ids, region, binsize))


def region_to_extent(h5, chrom_ids, region, binsize=None):
    return tuple(_region_to_extent(h5, chrom_ids, region, binsize))


class Cooler(object):
    """
    Convenient interface to a cooler HDF5 file.

    Parameters
    ----------
    store : str, h5py.File or h5py.Group
        File path or open handle to the root HDF5 group of a cooler.
    root : str, optional
        HDF5 path to root of cooler tree. If not provided, the root is assumed
        to be '/' if ``store`` is a file path or inferred if ``store`` is a 
        File or Group object.
    kwargs : optional
        Options to be passed to h5py.File() on every access. By default, the
        file is opened with the default driver and mode='r'.

    Notes
    -----
    * Metadata is accessible as a dictionary through the ``info`` property.

    * Data table range queries are provided through table selectors:
      ``chroms``, ``bins``, and ``pixels``, which return DataFrames
      or Series.

    * Matrix range queries are provided via a matrix selector, ``matrix``,
      which return ``scipy.sparse.coo_matrix`` arrays.

    If ``store`` is a file path, the file will be opened temporarily in
    read-only mode when performing operations. Such ``Cooler`` objects can be
    serialized and used multiprocessing, for example. See the following
    `discussion <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>`_
    on using h5py with multiprocessing safely.

    """
    def __init__(self, store, root=None, **kwargs):
        self.store = store
        self.kwargs = kwargs
        with open_hdf5(self.store, **self.kwargs) as h5:
            self.filename = h5.file.filename
            self.root = root if root is not None else h5.name
            grp = h5[self.root]

            _ct = chroms(grp)
            self._chromsizes = _ct.set_index('name')['length']
            self._chromids = dict(zip(_ct['name'], range(len(_ct))))
            self._info = info(grp)

    def _load_dset(self, path):
        with open_hdf5(self.store, **self.kwargs) as h5:
            grp = h5[self.root]
            return grp[path][:]

    def _load_attrs(self, path):
        with open_hdf5(self.store, **self.kwargs) as h5:
            grp = h5[self.root]
            return dict(grp[path].attrs)

    def open(self, mode='r', **kwargs):
        return h5py.File(self.store, mode, **kwargs)

    @property
    def chromsizes(self):
        """ Ordered mapping of reference sequences to their lengths in bp """
        return self._chromsizes

    @property
    def chromnames(self):
        """ List of reference sequence names """
        return list(self._chromsizes.index)

    def offset(self, region):
        """ Bin ID containing the left end of a genomic region

        Parameters
        ----------
        region : str or tuple
            Genomic range

        Returns
        -------
        int

        Examples
        --------
        >>> c.offset('chr3')  # doctest: +SKIP
        1311

        """
        with open_hdf5(self.store, **self.kwargs) as h5:
            grp = h5[self.root]
            return region_to_offset(
                grp, self._chromids,
                parse_region(region, self._chromsizes))

    def extent(self, region):
        """ Bin IDs containing the left and right ends of a genomic region

        Parameters
        ----------
        region : str or tuple
            Genomic range

        Returns
        -------
        2-tuple of ints

        Examples
        --------
        >>> c.extent('chr3')  # doctest: +SKIP
        (1311, 2131)

        """
        with open_hdf5(self.store, **self.kwargs) as h5:
            grp = h5[self.root]
            return region_to_extent(
                grp, self._chromids,
                parse_region(region, self._chromsizes))

    @property
    def info(self):
        """ File information and metadata

        Returns
        -------
        dict

        """
        with open_hdf5(self.store, **self.kwargs) as h5:
            grp = h5[self.root]
            return info(grp)

    @property
    def shape(self):
        return (self._info['nbins'],) * 2

    def chroms(self, **kwargs):
        """ Chromosome table selector

        Returns
        -------
        Table selector

        """
        def _slice(fields, lo, hi):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                return chroms(grp, lo, hi, fields, **kwargs)

        return RangeSelector1D(None, _slice, None, self._info['nchroms'])

    def bins(self, **kwargs):
        """ Bin table selector

        Returns
        -------
        Table selector

        """

        def _slice(fields, lo, hi):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                return bins(grp, lo, hi, fields, **kwargs)

        def _fetch(region):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                return region_to_extent(grp, self._chromids,
                                        parse_region(region, self._chromsizes))

        return RangeSelector1D(None, _slice, _fetch, self._info['nbins'])

    def pixels(self, join=False, **kwargs):
        """ Pixel table selector

        Parameters
        ----------
        join : bool, optional
            Whether to expand bin ID columns into chrom, start, and end columns.
            Default is ``False``.

        Returns
        -------
        Table selector

        """

        def _slice(fields, lo, hi):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                return pixels(grp, lo, hi, fields, join, **kwargs)

        def _fetch(region):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                i0, i1 = region_to_extent(
                    grp, self._chromids,
                    parse_region(region, self._chromsizes))
                lo = grp['indexes']['bin1_offset'][i0]
                hi = grp['indexes']['bin1_offset'][i1]
                return lo, hi

        return RangeSelector1D(None, _slice, _fetch, self._info['nnz'])

    def matrix(self, field=None, dense=False, balance=True, as_pixels=False, 
               join=False, ignore_index=True, max_chunk=500000000):
        """ Contact matrix selector

        Parameters
        ----------
        field : str, optional
            Which column of the pixel table to fill matrix selections with.
            By default, the 'count' column is used.
        balance : bool, optional
            Whether to apply pre-calculated matrix balancing weights to
            selections. Default is False.
        as_pixels : bool, optional
            Instead of a complete rectangular sparse matrix, return a DataFrame
            containing the corresponding rows from the pixel table.
            Default is False.
        join : bool, optional
            Whether to expand bin ID columns. False by default. Only applies if
            ``as_pixels`` is True.
        ignore_index : bool, optional
            If requesting pixels, don't populate the index column with the pixel
            IDs to improve performance. Default is True.

        Returns
        -------
        Matrix selector

        """

        def _slice(field, i0, i1, j0, j1):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                return matrix(grp, i0, i1, j0, j1, field, dense,
                    balance, as_pixels, join, ignore_index, max_chunk)

        def _fetch(region, region2=None):
            with open_hdf5(self.store, **self.kwargs) as h5:
                grp = h5[self.root]
                if region2 is None:
                    region2 = region
                region1 = parse_region(region, self._chromsizes)
                region2 = parse_region(region2, self._chromsizes)
                i0, i1 = region_to_extent(grp, self._chromids, region1)
                j0, j1 = region_to_extent(grp, self._chromids, region2)
                return i0, i1, j0, j1

        return RangeSelector2D(field, _slice, _fetch, (self._info['nbins'],) * 2)

    def __repr__(self):
        filename = os.path.basename(self.filename)
        return '<Cooler "{}":"{}">'.format(filename, self.root)


def info(h5):
    """
    File and user metadata dict.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.

    Returns
    -------
    dict

    """
    d = {}
    for k, v in h5.attrs.items():
        if isinstance(v, six.string_types):
            try:
                v = json.loads(v)
            except ValueError:
                pass
        d[k] = v
    return d


def chroms(h5, lo=0, hi=None, fields=None, **kwargs):
    """
    Table describing the chromosomes/scaffolds/contigs used.
    They appear in the same order they occur in the heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = (pandas.Index(['name', 'length'])
                        .append(pandas.Index(h5['chroms'].keys()))
                        .drop_duplicates())
    return get(h5['chroms'], lo, hi, fields, **kwargs)


def bins(h5, lo=0, hi=None, fields=None, **kwargs):
    """
    Table describing the genomic bins that make up the axes of the heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = (pandas.Index(['chrom', 'start', 'end'])
                        .append(pandas.Index(h5['bins'].keys()))
                        .drop_duplicates())
    return get(h5['bins'], lo, hi, fields, **kwargs)


def pixels(h5, lo=0, hi=None, fields=None, join=True, **kwargs):
    """
    Table describing the nonzero upper triangular pixels of the Hi-C contact
    heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.
    join : bool, optional
        Whether or not to expand bin ID columns to their full bin description
        (chrom, start, end). Default is True.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = (pandas.Index(['bin1_id', 'bin2_id', 'count'])
                        .append(pandas.Index(h5['pixels'].keys()))
                        .drop_duplicates())

    df = get(h5['pixels'], lo, hi, fields, **kwargs)

    if join:
        bins = get(h5['bins'], 0, None, ['chrom', 'start', 'end'])
        df = annotate(df, bins)

    return df


def annotate(pixels, bins, replace=True):
    """
    Add bin annotations to a selection of pixels by performing a "left join"
    from the bin IDs onto a table that describes properties of the genomic
    bins. New columns will be appended on the left of the output data frame.

    Parameters
    ----------
    pixels : DataFrame
        A data frame containing columns named ``bin1_id`` and/or ``bin2_id``.
        If columns ``bin1_id`` and ``bin2_id`` are both present in ``pixels``,
        the adjoined columns will be suffixed with '1' and '2' accordingly.
    bins : DataFrame or DataFrame view
        Data structure that contains a full description of the genomic bins of
        the contact matrix, where the index corresponds to bin IDs.
    replace : bool, optional
        Whether to remove the original ``bin1_id`` and ``bin2_id`` columns from
        the output. Default is True.

    Returns
    -------
    DataFrame

    """
    ncols = len(pixels.columns)

    if 'bin1_id' in pixels:
        if len(bins) > len(pixels):
            bin1 = pixels['bin1_id']
            lo = bin1.min()
            hi = bin1.max() + 1
            lo = 0 if np.isnan(lo) else lo
            hi = 0 if np.isnan(hi) else hi
            right = bins[lo:hi]
        else:
            right = bins[:]

        pixels = pixels.merge(
            right,
            how='left',
            left_on='bin1_id',
            right_index=True)

    if 'bin2_id' in pixels:
        if len(bins) > len(pixels):
            bin2 = pixels['bin2_id']
            lo = bin2.min()
            hi = bin2.max() + 1
            lo = 0 if np.isnan(lo) else lo
            hi = 0 if np.isnan(hi) else hi
            right = bins[lo:hi]
        else:
            right = bins[:]

        pixels = pixels.merge(
            right,
            how='left',
            left_on='bin2_id',
            right_index=True,
            suffixes=('1', '2'))

    # rearrange columns
    pixels = pixels[list(pixels.columns[ncols:]) + list(pixels.columns[:ncols])]

    # drop bin IDs
    if replace:
        cols_to_drop = [col for col in ('bin1_id', 'bin2_id') if col in pixels]
        pixels = pixels.drop(cols_to_drop, axis=1)

    return pixels


def matrix(h5, i0, i1, j0, j1, field=None, dense=False, balance=True,
           as_pixels=False, join=True, ignore_index=True, max_chunk=500000000):
    """
    Two-dimensional range query on the Hi-C contact heatmap.
    Returns either a rectangular sparse ``coo_matrix`` or a data frame of upper
    triangle pixels.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    i0, i1 : int, optional
        Bin range along the 0th (row) axis of the heatap.
    j0, j1 : int, optional
        Bin range along the 1st (col) axis of the heatap.
    field : str, optional
        Which column of the pixel table to fill the matrix with. By default,
        the 'count' column is used.
    balance : bool, optional
        Whether to apply pre-calculated matrix balancing weights to the
        selection. Default is False.
    as_pixels: bool, optional
        Return a DataFrame of the corresponding rows from the pixel table
        instead of a rectangular sparse matrix. False by default.
    join : bool, optional
        If requesting pixels, specifies whether to expand the bin ID columns
        into (chrom, start, end). Has no effect when requesting a rectangular
        matrix. Default is True.
    ignore_index : bool, optional
        If requesting pixels, don't populate the index column with the pixel
        IDs to improve performance. Default is True.

    Returns
    -------
    coo_matrix

    Notes
    -----
    Use the ``toarray()`` method to convert to a sparse matrix to a dense
    NumPy array.

    """
    if field is None:
        field = 'count'

    triu_reader = TriuReader(h5, field, max_chunk)

    if balance and 'weight' not in h5['bins']:
        raise ValueError(
            "No column 'bins/weight' found. Use ``cooler.ice`` to "
            "calculate balancing weights or set balance=False.")

    if as_pixels:
        index = None if ignore_index else triu_reader.index_col(i0, i1, j0, j1)
        i, j, v = triu_reader.query(i0, i1, j0, j1)

        cols = ['bin1_id', 'bin2_id', field]
        df = pandas.DataFrame(dict(zip(cols, [i, j, v])),
                              columns=cols, index=index)

        if balance:
            weights = get(h5['bins'], min(i0, j0), max(i1, j1), ['weight'])
            df2 = annotate(df, weights)
            df['balanced'] = df2['weight1'] * df2['weight2'] * df2[field]

        if join:
            bins = get(h5['bins'], min(i0, j0), max(i1, j1), ['chrom', 'start', 'end'])
            df = annotate(df, bins)

        return df

    elif dense:
        i, j, v = query_rect(triu_reader.query, i0, i1, j0, j1)
        arr = coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0)).toarray()

        if balance:
            weights = h5['bins']['weight']
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            arr = arr * np.outer(bias1, bias2)

        return arr

    else:
        i, j, v = query_rect(triu_reader.query, i0, i1, j0, j1)
        mat = coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0))

        if balance:
            weights = h5['bins']['weight']
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            mat.data = bias1[mat.row] * bias2[mat.col] * mat.data

        return mat
