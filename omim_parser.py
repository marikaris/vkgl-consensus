class OmimParser:
    def __init__(self, filename):
        self.codes = self.parse(filename)

    def parse(self, filename):
        omim_genes = {}
        for line in open(filename):
            vals = line.split('\t')
            if vals[1].strip('\n') in omim_genes:
                omim_genes[vals[1].strip('\n')]+= ','+vals[0]
            else:
                omim_genes[vals[1].strip('\n')] = vals[0]
        return omim_genes