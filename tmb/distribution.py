from tmb.tmbresult import TMBResult
from tmb.bamreader import BamReader
from tmb.config import Config
from scipy import stats
from scipy import spatial
import logging
import math


def get_tmb_from_files(tumor_path, normal_path, exon):
    logging.info(f'Collect distributions from {exon.chr}:{exon.start}-{exon.end}')
    tumor_distribution, normal_distribution = collect_distributions(tumor_path, normal_path, exon)

    logging.info(f'Hypotesis test for {exon.chr}:{exon.start}-{exon.end}')
    somatic_sites = hypotesis_test(tumor_distribution, normal_distribution)

    logging.info(f'Calculates TMB for {exon.chr}:{exon.start}-{exon.end}')
    tmb = calculate_tmb_from_sites(exon, somatic_sites)

    logging.info(f'Finished calculation of TMB for {exon.chr}:{exon.start}-{exon.end}')
    return tmb


def calculate_tmb_from_sites(exon, somatic_sites):
    tmb_value = round((len(somatic_sites) / exon.length) * 1000000, 2)
    tmb = TMBResult(exon, tmb_value, somatic_sites)
    return tmb


def collect_distributions(tumor_path, normal_path, exon):
    tumor = BamReader(tumor_path).file
    normal = BamReader(normal_path).file

    tumor_distribution = fill_maps_reads(tumor, exon)
    tumor.close()
    # print(tumor_distribution)
    normal_distribution = fill_maps_reads(normal, exon)
    normal.close()
    # print(normal_distribution)
    return tumor_distribution, normal_distribution


def fill_maps_pileup(bam, exon):
    positions_to_acgt = {}
    for i in range(exon.start, exon.end +1):
        acgt_to_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "DEL": 0, "INS": 0}
        positions_to_acgt[i] = acgt_to_counts

    # Pile can get reads that will cover this position in column-like form
    pile = bam.pileup(exon.chr, exon.start, exon.end,
                      stepper="nofilter",
                      ignore_overlaps=False,
                      flag_filter=0,
                      ignore_orphans=False)

    for pileupcolumn in pile:
        position = pileupcolumn.pos
        if position < exon.start - 1 or position >= exon.end:
            continue

        # print(f'position: {position}')
        good_reads = list(filter(good_pileupread, pileupcolumn.pileups))
        for pileupread in good_reads:
            if pileupread.indel > 0:
                insertion(exon, pileupread, position, positions_to_acgt)
            if pileupread.indel < 0 or pileupread.query_position is None:
                deletion(exon, pileupread, position, positions_to_acgt)
            else:
                single_nucleotide(pileupread, position, positions_to_acgt)
    return positions_to_acgt


def fill_maps_reads(bam, exon):
    positions_to_acgt = {}
    for i in range(exon.start, exon.end + 1):
        acgt_to_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "DEL": 0, "INS": 0}
        positions_to_acgt[i] = acgt_to_counts

    fetch = bam.fetch(exon.chr, exon.start, exon.end)

    for read in fetch:
        start = read.reference_start
        end = read.reference_end
        position = 0
        opers = 0
        if good_read(read) and read.cigartuples:
            for oper in read.cigartuples:
                if opers + start > exon.end:
                    break

                # M
                if oper[0] == 0:
                    opers += oper[1]
                    shift, next = check_oper(opers, start, exon, read)
                    if next:
                        continue
                    i = 0
                    while i < oper[1] - shift and start + position + shift + i < exon.end:
                        base_quality = read.query_qualities[position + shift + i]
                        if base_quality >= Config.baseq:
                            base = read.query_sequence[position + shift + i]
                            positions_to_acgt[start + position + shift + i + 1][base] += 1
                        i += 1
                    position += oper[1]

                # S, H or N
                if oper[0] == 3 or oper[0] == 4 or oper[0] == 5:
                    opers += oper[1]
                    position += oper[1]
                    shift, next = check_oper(opers, start, exon, read)
                    if next:
                        continue
                # I
                if oper[0] == 1:
                    opers += oper[1]
                    shift, next = check_oper(opers, start, exon, read)
                    if next:
                        continue
                    i = 0
                    while i < oper[1]:
                        positions_to_acgt[start + position + i + 1]["INS"] += 1
                        i += 1
                        #position += 1
                # D
                if oper[0] == 2:
                    opers += oper[1]
                    shift, next = check_oper(opers, start, exon, read)
                    if next:
                        continue
                    i = 0
                    while i < oper[1]:
                        positions_to_acgt[start + position + i + 1]["DEL"] += 1
                        i += 1
                        #position += 1

    return positions_to_acgt


def check_oper(opers, start, exon, read):
    shift = 0
    if opers + start >= exon.start:
        shift = read.query_length - (opers + start - exon.start) - 1

    if opers + start < exon.start:
        return shift, True
    return shift, False

def single_nucleotide(pileupread, position, positions_to_acgt):
    base_quality = pileupread.alignment.query_qualities[pileupread.query_position]
    if base_quality >= Config.baseq:
        base = pileupread.alignment.query_sequence[pileupread.query_position]
        positions_to_acgt[position + 1][base] += 1


def deletion(exon, pileupread, position, positions_to_acgt):
    # print("DEL", pileupread.indel)
    bases = 0
    while position + 1 - bases + 1 <= exon.end and bases > pileupread.indel:
        base = 'DEL'
        positions_to_acgt[position + 1 - bases + 1][base] += 1
        bases = bases - 1


def insertion(exon, pileupread, position, positions_to_acgt):
    # print("INS", pileupread.indel)
    bases = 1
    while position + bases + 1 <= exon.end and bases <= pileupread.indel:
        base = 'INS'
        positions_to_acgt[position + bases + 1][base] += 1
        bases = bases + 1


def hypotesis_test(tumor, normal):
    somatic_sites = []
    for position in tumor:
        tumor_bases = tumor[position]
        normal_bases = normal[position]

        tumor_counts = list(tumor_bases.values())
        normal_counts = list(normal_bases.values())

        if not bad_counts(normal_counts, tumor_counts):
            cosine = dist_by_cosine(normal_counts, tumor_counts)
            #dist_by_percentage(normal_counts, tumor_counts, position, somatic_sites)

            if cosine > 0.005:
                somatic_sites.append(position)
    return somatic_sites


def bad_counts(normal_counts, tumor_counts):
    return sum(tumor_counts) == 0 \
           or sum(normal_counts) == 0 \
           or max(tumor_counts) < Config.mincov \
           or max(normal_counts) < Config.mincov


# Not used
def adjust_zero_counts(normal_counts, tumor_counts):
    # To adjust zero coverages so we can use chi-test
    tumor_total_coverage = sum(tumor_counts)
    normal_total_coverage = sum(normal_counts)
    tumor_percentage = list([x * 100 / tumor_total_coverage for x in tumor_counts])
    normal_percentage = list([x * 100 / normal_total_coverage for x in normal_counts])
    return normal_percentage, tumor_percentage


def dist_by_pearsonr_test(normal, tumor, position, somatic_sites):
    pearsonr = stats.pearsonr(tumor, normal)
    if pearsonr[1] > 0.005:
        # Sam is 0-based, extend position:
        somatic_sites.append(position)


def dist_by_eucleadin(normal, tumor):
    sum = 0
    for i in range(len(tumor)):
        dist = (normal[i] - tumor[i])**2
        sum += dist
    return


def dist_by_cosine(normal, tumor):
    # print(normal)
    # print(tumor)
    cosine = spatial.distance.cosine(normal, tumor)

    return cosine


def dist_by_percentage(normal_counts, tumor_counts, position, somatic_sites):
    for nm, tm in zip(normal_counts, tumor_counts):
        distance = abs(nm/sum(normal_counts) - tm/sum(tumor_counts))
        if distance > Config.freq and tm > Config.mincov and nm > Config.mincov:
            # Sam is 0-based, extend position:
            somatic_sites.append(position)
            break


def good_pileupread(read):
    """good_read(read) --> Check if read fits criteria for quality

    :param read: SAM file record
    :return: true if read is fit criteria, false if not
    """
    if read.alignment.mapping_quality < Config.mapq:
        return False
    # TODO: add other checks
    return True


def good_read(read):
    """good_read(read) --> Check if read fits criteria for quality

    :param read: SAM file record
    :return: true if read is fit criteria, false if not
    """
    if read.mapping_quality < Config.mapq:
        return False
    # TODO: add other checks
    return True