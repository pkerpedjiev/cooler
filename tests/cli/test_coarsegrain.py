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
from cooler.cli.coarsegrain import coarsegrain


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
        coarsegrain, [
            op.join(testdir, 'data', 'dec2_20_pluslig_1pGene_grch38_UBR4_D_1nt.pairwise.sorted.cool'),
            '--output-file', multires_path,
            '--no-balance',
        ]
    )

    sys.stdout.write(result.output)

    # this file should have 6 zoom levels

    #c = cooler.Cooler(h5py.File(multires_path)['5'])
    #print('pixels (5):', len(c.pixels()[:].index))

    # should get a ValueError because the chromosome names in the pixels dont' match
    # the stored chromosome names
    assert(result.exit_code == -1)

@with_setup(teardown=partial(teardown_func, multires_path))
def test_coarsegrain():
    runner = CliRunner()
    result = runner.invoke(
        coarsegrain, [
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
            '--output-file', multires_path
        ]
    )

    sys.stdout.write(result.output)

    # this file should have 6 zoom levels

    #c = cooler.Cooler(h5py.File(multires_path)['5'])
    #print('pixels (5):', len(c.pixels()[:].index))

    # should get a ValueError because the chromosome names in the pixels dont' match
    # the stored chromosome names
    assert(result.exit_code == 0)

    import cooler.contrib.higlass as cch

    with h5py.File(multires_path) as f:
        d = cch.get_data(f, 0, 0, 3000000000, 0, 3000000000)
        print("d:", d.iloc[0]['balanced'])
        assert(d.iloc[0]['balanced'] -  0.8179 < 0.001)

        d = cch.get_data(f, 0, 1000000000, 3000000000, 1000000000, 3000000000)
        print("d:", d.iloc[0]['balanced'])
        assert(d.iloc[0]['balanced'] -  0.8179 > 0.001)

    runner = CliRunner()
    result = runner.invoke(
        coarsegrain, [
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
            '--output-file', multires_path,
            '--offset-x', 1000000000,
            '--offset-y', 2000000000,
        ]
    )

    sys.stdout.write(result.output)

    # this file should have 6 zoom levels

    #c = cooler.Cooler(h5py.File(multires_path)['5'])
    #print('pixels (5):', len(c.pixels()[:].index))

    # should get a ValueError because the chromosome names in the pixels dont' match
    # the stored chromosome names
    assert(result.exit_code == 0)
