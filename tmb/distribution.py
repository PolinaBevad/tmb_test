from tmb.tmbresult import TMBResult
from tmb.bamreader import BamReader
from tmb.config import Config


def calculate_tmb(tumor_path, normal_path, exon):
    tumor_distribution, normal_distribution = collect_distributions(tumor_path, normal_path, exon)
    tmb, somatic_sites = hypotesis_test(tumor_distribution, normal_distribution, exon)
    tmb = TMBResult(exon, tmb, somatic_sites)
    return tmb


def collect_distributions(tumor_path, normal_path, exon):
    tumor = BamReader(tumor_path).file
    normal = BamReader(normal_path).file

    tumor_distribution = fill_maps(tumor, exon)
    # print(tumor_distribution)
    normal_distribution = fill_maps(normal, exon)
    # print(normal_distribution)
    return tumor_distribution, normal_distribution


def fill_maps(bam, exon):
    positions_to_acgt = {}
    for i in range(exon.start, exon.end):
        positions_to_acgt[i] = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "DEL": 0}

    # Pile can get reads that will cover this position in column-like form
    pile = bam.pileup(exon.chr, exon.start, exon.end)
    for pileupcolumn in pile:
        position = pileupcolumn.pos
        if position < exon.start or position >= exon.end:
            continue

        acgt_to_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "DEL": 0}
        for pileupread in pileupcolumn.pileups:
            if bad_read(pileupread):
                continue
            if pileupread.query_position is None:
                base = 'DEL'
            else:
                base_quality = pileupread.alignment.query_qualities[pileupread.query_position]
                if base_quality < Config.baseq:
                    break
                base = pileupread.alignment.query_sequence[pileupread.query_position]
            acgt_to_counts[base] += 1
        positions_to_acgt[position] = acgt_to_counts
    return positions_to_acgt


def hypotesis_test(tumor, normal, exon):
    exone_len = exon.end - exon.start
    somatic_sites = []
    for position in tumor:
        tumor_counts = list(tumor[position].values())
        normal_counts = list(normal[position].values())

        tumor_total_coverage = sum(tumor_counts)
        normal_total_coverage = sum(normal_counts)

        if tumor_total_coverage == 0 or normal_total_coverage == 0:
            continue

        tumor_percentage = list([x / tumor_total_coverage for x in tumor_counts])
        normal_percentage = list([x / normal_total_coverage for x in normal_counts])

        for i in range(len(tumor_counts)):
            a = abs(tumor_percentage[i] - normal_percentage[i])
            if a > Config.freq and tumor_counts[i] > Config.mincov and normal_counts[i] > Config.mincov:
                # Sam is 0-based, extend position:
                somatic_sites.append(position + 1)
                break

    tmb = round((len(somatic_sites) / exone_len) * 1000000, 2)
    return tmb, somatic_sites


# Read doesn't fit criteria for quality
def bad_read(read):
    if read.alignment.mapping_quality < Config.mapq:
        return True

    return False

