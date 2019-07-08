#!/usr/bin/env python

import argparse
from tmb.config import Config
from tmb.bedreader import BedReader
from tmb.distribution import get_tmb_from_files
from multiprocessing import Pool
import logging


def main():
    config = parse_config()
    set_logging(config)
    exons = BedReader(config.bed).exonList

    pool = Pool(processes=Config.th)
    results = [pool.apply_async(get_tmb_from_files, args=(config.tumor, config.normal, exon)) for exon in exons]

    if Config.header:
        print_header()

    for output in results:
        tmb = output.get()
        print(tmb)


def parse_config():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required")
    optional = parser.add_argument_group("optional")
    optional.add_argument('--freq', type=float,
                        help="Minimum difference of frequencies for base to consider site as somatic. "
                             "Default: 5%% (0.05)")
    optional.add_argument('--mapq', type=float, help="Minimum read mapping quality. Default: 10.0")
    optional.add_argument('--baseq', type=int, help="Minimum base quiality. Default: 25")
    optional.add_argument('--mincov', type=int,
                        help="Minimum coverage of base for nucleotide to be considered as somatic. "
                             "Default: 2")
    optional.add_argument('--th', type=int, help="Number of threads for multiprocessing mode. Default: 1")
    optional.add_argument('--pvalue', type=int, help="Pvalue threshold for somatic sites. Default: 5*10E20")
    optional.add_argument('--log', type=str, help="The level of logging from CRITICAL, ERROR, WARNING, INFO, DEBUG. ")
    optional.add_argument('--header', type=str, help="Print header to the output")
    required.add_argument('--normal', required=True, type=str, help="Path to normal SAM/BAM/CRAM file.")
    required.add_argument('--tumor', required=True, type=str, help="Path to tumor SAM/BAM/CRAM file.")
    required.add_argument('--bed', required=True, type=str, help="Path to BED file.")

    args = parser.parse_args()
    return Config(args)


def print_header():
    print("\t".join(['Chr', 'Start', 'End', 'Gene', 'SomaticSites', 'TMB']))


def set_logging(config):
    log_level = getattr(logging, Config.log_level.upper(), None)
    if not isinstance(log_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(filename='tmb.log', level=log_level, filemode='w',
                        format='%(asctime)s %(levelname)s:%(message)s ')
    logging.info("Starting TMB.")
    logging.info("Retrieving exons from %s.", config.bed)


if __name__ == "__main__":
    main()
