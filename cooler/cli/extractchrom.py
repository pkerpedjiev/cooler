# -*- coding: utf-8 -*-
from __future__ import division, print_function
from six.moves import map
import multiprocess as mp
import tempfile
import argparse
import logging
import sys

import numpy as np
import h5py

import cooler

from ..io import CoolerAggregator, create
from ..ice import iterative_correction
from ..util import binnify, get_chromsizes
from ..tools import lock
from .. import api

import click
from . import cli, logger

@cli.command()
@click.argument(
    'cooler_file',
    metavar="COOLER_PATH")
@click.argument(
    'chromosome',
    metavar='CHROMOSOME')
@click.argument(
    'output_file',
    metavar="OUTPUT_FILE")
def extractchrom(cooler_file, chromosome, output_file):
    """
    Extract a chromosome from a cooler file

    cooler_file: The input file
    chromosome: The name of the chromosome to keep
    output_file: The name of the output file
    """
    infile = cooler_file
    chunksize = 1000000
    output_bin1_id = []
    output_bin2_id = []
    output_count = []

    with h5py.File(cooler_file, 'r') as src:
        with h5py.File(output_file, 'w') as dest:
            src.copy('/bins', dest, 'bins') 
            src.copy('/chroms', dest, 'chroms') 
            src.copy('/indexes', dest, 'indexes')

            for key,value in src.attrs.items():
                dest.attrs[key] = value

            c = cooler.Cooler(src)
            pixels_length = len(src['pixels']['bin1_id'])
            print("pixels_length:", pixels_length)

            # the total number of pixels seen so far
            count = 0
            # the total number of pixels in the chosen chromosome
            chrom_count = 0

            # do an initial pass to count how many times this chromosome occurs
            # in the dataset
            while count < pixels_length:
                chunk = min(pixels_length - count, chunksize)

                pixels = c.pixels()[count:count+chunk]
                joined_pixels = c.pixels(join=True)[count:count+chunk]

                chrom_count += sum(((joined_pixels['chrom1'] == chromosome) & (joined_pixels['chrom2'] == chromosome)).values)
                count += chunk
                print("First pass (out of two), pixels read:", count)

            count = 0

            grp = dest.create_group('pixels')

            # the current position of the kept bin ids
            current_kept_pos = 0

            pixels_bin1_ids = grp.create_dataset('bin1_id', (chrom_count,), dtype='i8')
            pixels_bin2_ids = grp.create_dataset('bin2_id', (chrom_count,), dtype='i8')
            pixels_count = grp.create_dataset('count', (chrom_count,), dtype='i8')

            while count < pixels_length:
                chunk = min(pixels_length - count, chunksize)

                pixels = c.pixels()[count:count+chunk]
                joined_pixels = c.pixels(join=True)[count:count+chunk]
                to_keep_positions = (joined_pixels['chrom1'] == chromosome) & (joined_pixels['chrom2'] == chromosome)
                to_keep_pixels = pixels[to_keep_positions]

                num_to_keep = len(to_keep_pixels)
                print('num_to_keep', num_to_keep)

                if num_to_keep > 0:
                    pixels_bin1_ids[current_kept_pos:current_kept_pos + num_to_keep] = to_keep_pixels['bin1_id'].values
                    pixels_bin2_ids[current_kept_pos:current_kept_pos + num_to_keep] = to_keep_pixels['bin2_id'].values
                    pixels_count[current_kept_pos:current_kept_pos + num_to_keep] = to_keep_pixels['count'].values

                    current_kept_pos += num_to_keep

                count += chunk
                print("Second pass (out of two), pixels read:", count)
