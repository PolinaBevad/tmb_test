#!/usr/bin/env python

import argparse
from tmb.config import Config
from tmb.bedreader import BedReader
from tmb.distribution import calculate_tmb
from multiprocessing import Pool


#TODO: add logger
def main():
    config = parse_config()
    exons = BedReader(config.bed).exonList

    pool = Pool(processes=Config.th)
    results = [pool.apply_async(calculate_tmb, args=(config.tumor, config.normal, exon)) for exon in exons]

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
    required.add_argument('--normal', required=True, type=str, help="Path to normal SAM/BAM/CRAM file.")
    required.add_argument('--tumor', required=True, type=str, help="Path to tumor SAM/BAM/CRAM file.")
    required.add_argument('--bed', required=True, type=str, help="Path to BED file.")

    args = parser.parse_args()
    return Config(args)


def print_header():
    print("\t".join(['Chr', 'Start', 'End', 'Gene', 'SomaticSites', 'TMB']))


if __name__ == "__main__":
    main()
