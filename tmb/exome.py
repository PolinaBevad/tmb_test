class Exome:
    def __init__(self, chr, start, end, gene):
        self.chr = chr
        self.start = start
        self.end = end
        self.gene = gene
        self.length = end - start

    def __str__(self):
        exon_description = ', '.join(['{key}={value}'.format(key=key, value=self.__dict__.get(key))
                                      for key in self.__dict__])
        return '\n' + exon_description




