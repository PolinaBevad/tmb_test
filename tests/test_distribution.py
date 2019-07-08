from numpy.testing import assert_almost_equal

from tmb.bedreader import BedReader
from tmb.bamreader import BamReader
from tmb.exome import Exome
from tmb.distribution import get_tmb_from_files
from tmb.distribution import calculate_tmb_from_sites
from tmb.distribution import fill_maps_reads
from tmb.distribution import fill_maps_pileup
from tmb.distribution import hypotesis_test

import pytest


def test_distribution_by_pileup():
    normal = '../data/chr5_665281/normal_chr5_665281.bam'
    tumor = '../data/chr5_665281/tumour_chr5_665281.bam'
    exons = BedReader('../data/chr5_665281/dist.bed').exonList

    for exon in exons:
        tmb = get_tmb_from_files(tumor, normal, exon)
        assert tmb.somatic_sites == [665281, 665308, 665336]
        assert_almost_equal(tmb.tmb, 50847.46)


def test_panel_indexes():
    tumor = '../data/panel_az_600/Dev_731_GTL_16_5_Pool1-ready.bam'
    normal = '../data/panel_az_600/Dev_731_NA12878a_Pool1-ready.bam'
    exons = BedReader('../data/panel_az_600/panel_az600_chr7_MET.bed').exonList

    for exon in exons:
        tmb = get_tmb_from_files(tumor, normal, exon)
        print(tmb)


def test_cosine_and_other_dist():
    exon = Exome('chr1', 0, 2, "")
    normal = {0: {'T': 77, 'DEL': 19}}
    tumor = {0: {'T': 23, 'DEL': 10}}
    somatic_sites = hypotesis_test(normal, tumor)
    tmb = calculate_tmb_from_sites(exon, somatic_sites)

    assert len(tmb.somatic_sites) == 1
    assert tmb.somatic_sites == [0]
    assert tmb.tmb == 500000


def test_fill_maps_reads():
    normal_file = '../data/chr5_665281/normal_chr5_665281.bam'
    tumor_file = '../data/chr5_665281/tumour_chr5_665281.bam'
    exons = BedReader('../data/chr5_665281/dist.bed').exonList

    tumor = BamReader(tumor_file).file
    normal = BamReader(normal_file).file

    for exon in exons:
        mapNormalR = fill_maps_reads(normal, exon)
        mapNormalP = fill_maps_pileup(normal, exon)
        mapTumorR = fill_maps_reads(tumor, exon)
        mapTumorP = fill_maps_pileup(tumor, exon)

        print("normal reads:", mapNormalR)
        print("normal pileup:", mapNormalP)

        print("  tumor reads:", mapTumorR)
        print(" tumor pileup:", mapTumorP)


def test_string():
    str = 'AAGTGCGAGGTGACAGAGTCCCTGCTCTGGGGACACCGGCAGTGCCGCGAGGCATGGAGCGGGGGCTACTTGCCCCCTTGCACCTTCTGGTCGTAGGTGCCTGCGTGCTTGTAGCCGGACACATAGCCAGACTCGTCCACCAGATCCACG'
    print(str[128])
