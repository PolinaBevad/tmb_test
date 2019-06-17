from tmb.bamreader import BamReader
from tmb.bedreader import BedReader

import pytest


def test_bam_test_read():
    path = '../data/bed_bam_test/test.bam'
    bam = BamReader(path).file
    bam_iter = bam.fetch('1', 28234090, 28234093)
    list1 = [x for x in bam_iter]
    assert len(list1) == 4


def test_bed_test_read_correct():
    path = '../data/bed_bam_test/test1.bed'
    exons = BedReader(path).exonList
    assert len(exons) == 4

    path = '../data/bed_bam_test/test2.bed'
    exons = BedReader(path).exonList
    assert len(exons) == 4


def test_bed_test_read_incorrect():
    path = '../data/bed_bam_test/test3.bed'
    with pytest.raises(IndexError):
        BedReader(path)

