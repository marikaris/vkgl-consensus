import pprint
class OmimParser:
    def __init__(self, filename):
        self.codes = self.parse(filename)

    def parse(self, filename):
        omim_genes = {}
        for line in open(filename):
            line = line.split('\t')
            if len(line[1]) > 0 and len(line[2]) > 2:
                if line[2].strip('\n') in omim_genes:
                    if line[1] not in omim_genes[line[2].strip('\n')]:
                        omim_genes[line[2].strip('\n')].append(line[1])
                else:
                    omim_genes[line[2].strip('\n')] = [line[1]]
        return omim_genes

if __name__ == '__main__':
    pprint.pprint(OmimParser('omim.txt').codes)