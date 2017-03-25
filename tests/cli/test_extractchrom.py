# -*- coding: utf-8 -*-
from functools import partial
import os.path as op
import subprocess
import tempfile
import filecmp
import sys
import os
import cooler
import logging

import numpy as np
import h5py

import nose
from nose.tools import with_setup, set_trace
from click.testing import CliRunner

from cooler.cli.cload import cload, tabix as cload_tabix
from cooler.cli.extractchrom import extractchrom


if sys.version_info[0] == 3 and sys.version_info[1] == 3:
    raise nose.SkipTest


testdir = op.realpath(op.join(op.dirname(__file__), op.pardir))
tmp = tempfile.gettempdir()

multires_path = op.join(tmp, 'test.multires.cool')


def teardown_func(*filepaths):
    for fp in filepaths:
        try:
            os.remove(fp)
        except OSError:
            pass


@with_setup(teardown=partial(teardown_func, multires_path))
def test_csort():
    runner = CliRunner()
    result = runner.invoke(
        extractchrom, [
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
            'chrX',
            multires_path
        ]
    )

    sys.stdout.write(result.output)
    sys.stdout.write(str(result.exception))

    with h5py.File(multires_path) as f:
        bins = f['bins']['chrom']
        print("bins:", bins)

    # this file should have 6 zoom levels

    #c = cooler.Cooler(h5py.File(multires_path)['5'])
    #print('pixels (5):', len(c.pixels()[:].index))

    # should get a ValueError because the chromosome names in the pixels dont' match
    # the stored chromosome names
    assert(result.exit_code == 0)
