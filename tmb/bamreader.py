import pysam


class BamReader:
    def __init__(self, path):
        extension = path.strip().split('.')[-1]
        if extension == 'bam':
            self.file = self.read_bam(path)
        if extension == 'cram':
            self.file = self.read_cram(path)
        if extension == 'cram':
            self.file = self.read_sam(path)

    def read_bam(self, path):
        file = pysam.AlignmentFile(path, "rb")
        return file

    def read_sam(self, path):
        file = pysam.AlignmentFile(path, "r")
        return file

    def read_cram(self, path):
        file = pysam.AlignmentFile(path, "rc")
        return file
