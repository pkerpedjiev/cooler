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
def test_extractchrom():
    runner = CliRunner()
    result = runner.invoke(
        extractchrom, [
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
            'chr14',
            multires_path
        ],
        catch_exceptions=False
    )

    #sys.stdout.write(result.output)
    #sys.stdout.write(str(result.exception))
    c = cooler.Cooler(multires_path)

    assert(len(c.pixels()) > 0)
    assert(result.exit_code == 0)

@with_setup(teardown=partial(teardown_func, multires_path, multires_path + ".multires"))
def test_extract_and_aggregate():
    runner = CliRunner()
    result = runner.invoke(
        extractchrom, [
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
            'chrX',
            multires_path
        ],
        catch_exceptions=True
    )

    #sys.stdout.write(result.output)
    #sys.stdout.write(str(result.exception))

    with h5py.File(multires_path) as f:
        bins = f['bins']['chrom']

    assert(result.exit_code == 0)

    result = runner.invoke(
            coarsegrain, [
                multires_path,
                '--output-file', multires_path + '.multires'],
            catch_exceptions=True)

    #sys.stdout.write(result.output)
    #sys.stdout.write(str(result.exception))

    assert(result.exit_code == 0)
