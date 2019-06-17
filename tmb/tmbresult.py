class TMBResult:
    def __init__(self, exon, tmb, somatic_sites):
        self.exon = exon
        self.tmb = tmb
        self.somatic_sites = somatic_sites

    def __str__(self):
        return "\t".join([str(self.exon.chr), str(self.exon.start), str(self.exon.end), self.exon.gene,
                         str(len(self.somatic_sites)), str(self.somatic_sites), str(self.tmb)])