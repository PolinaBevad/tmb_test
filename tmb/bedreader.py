from tmb.exome import Exome


class BedReader:

    def __init__(self, path):
        extension = path.split('.')[-1]
        if extension == 'bed':
            self.exonList = self.read_bed(path)

    def read_bed(self, path):
        exonlist = []

        file = open(path, "r")
        for x in file:
            line = x.split()
            if len(line) == 4:
                exome = Exome(str(line[0]), int(line[1]), int(line[2]), str(line[3]))
            else:
                exome = Exome(str(line[0]), int(line[1]), int(line[2]), '.')
            exonlist.append(exome)
        return exonlist
