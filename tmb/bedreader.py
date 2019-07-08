from tmb.exome import Exome


class BedReader:

    def __init__(self, path):
        extension = path.split('.')[-1]
        if extension == 'bed':
            self.exonList = self.read_bed(path)

    def read_bed(self, path):
        exonlist = []

        with open(path, "r") as file:
            for line in file:
                if not line == '\n':
                    line_cols = line.split()
                    if len(line_cols) == 4:
                        exome = Exome(str(line_cols[0]), int(line_cols[1]), int(line_cols[2]), str(line_cols[3]))
                    else:
                        exome = Exome(str(line_cols[0]), int(line_cols[1]), int(line_cols[2]), '.')
                    exonlist.append(exome)
        return exonlist
