from numpy.testing import assert_almost_equal

from tmb.bedreader import BedReader
from tmb.distribution import calculate_tmb

import pytest


def test_distribution():
    normal = '../data/chr5_665281/normal_chr5_665281.bam'
    tumor = '../data/chr5_665281/tumour_chr5_665281.bam'
    exons = BedReader('../data/chr5_665281/dist.bed').exonList

    for exon in exons:
        tmb = calculate_tmb(tumor, normal, exon)
        assert tmb.somatic_sites == [665281, 665308, 665336]
        assert_almost_equal(tmb.tmb, 50847.46)


def test_panel_indexes():
    tumor = '../data/panel_az_600/Dev_731_GTL_16_5_Pool1-ready.bam'
    normal = '../data/panel_az_600/Dev_731_NA12878a_Pool1-ready.bam'
    exons = BedReader('../data/panel_az_600/panel_az600_chr7_MET.bed').exonList

    for exon in exons:
        tmb = calculate_tmb(tumor, normal, exon)
