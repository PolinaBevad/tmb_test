# Configuration parameters for thresholds of quality etc.


class Config:
    freq = 0.05
    mapq = 10.0
    baseq = 25
    mincov = 2
    th = 1
    pvalue = 10E-60
    log_level = "WARNING"
    header = False

    @staticmethod
    def set_config(args):
        if args.freq:
            Config.freq = args.freq
        if args.mapq:
            Config.mapq = args.mapq
        if args.baseq:
            Config.baseq = args.baseq
        if args.mincov:
            Config.mincov = args.mincov
        if args.th:
            Config.th = args.th
        if args.pvalue:
            Config.pvalue = args.pvalue
        if args.log:
            Config.log_level = args.log
        if args.header:
            Config.header = True

    def __init__(self, args):
        self.normal = args.normal
        self.tumor = args.tumor
        self.bed = args.bed
        self.set_config(args)
