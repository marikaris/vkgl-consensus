class MolgenisConfigParser():
    def __init__(self, file):
        self.config = self.parse(file)

    def parse(self, file):
        config = {}
        for line in open(file):
            values = line.split('=')
            config[values[0]] = values[1]
        return config