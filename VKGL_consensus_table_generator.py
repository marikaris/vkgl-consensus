import molgenis
from Molgenis_config_parser import MolgenisConfigParser
import pprint

class ConsensusTableGenerator():
    def __init__(self, labs, session):
        self.labs = labs
        self.session = session
        self.old_diseases = {}
        self.old_comments = {}
        self.clear_tables()
        print('Consensus and consensus comments cleared')
        self.lab_data = self.process_data()
        table = self.calculate_consensus()
        self.upload_consensus(table)

    def process_data(self):
        consensus = {}
        for lab in self.labs:
            lab_data = self.session.get('VKGL_' + lab, num=10000)
            for variant in lab_data:
                variantId = variant['id'].replace(lab + '_', '')
                if variantId not in consensus:
                    protein = ['' if 'protein' not in variant else variant['protein']][0]
                    disease = ['' if 'consensus_' + variantId not in self.old_diseases else self.old_diseases[
                        'consensus_' + variantId]][0]
                    consensus[variantId] = {lab + '_classification': variant['id'], 'counter': {'b': 0, 'p': 0, 'v': 0},
                                            'REF': variant['REF'], 'ALT': variant['ALT'], 'gene': variant['gene'],
                                            'cDNA': variant['cDNA'], 'protein': protein,
                                            'chromosome': str(variant['chromosome']), 'stop': str(variant['stop']),
                                            'POS': str(variant['POS']), 'id': 'consensus_' + variantId,
                                            'comments': 'consensus_' + variantId,
                                            lab.lower(): variant['classification']}
                    if disease != '':
                        consensus[variantId]['disease'] = str(disease)
                else:
                    consensus[variantId][lab + '_classification'] = variant['id']
                    consensus[variantId][lab.lower()] = variant['classification']
                if variant['classification'] == 'Benign' or variant['classification'] == 'Likely benign':
                    consensus[variantId]['counter']['b'] += 1
                elif variant['classification'] == 'Pathogenic' or variant['classification'] == 'Likely pathogenic':
                    consensus[variantId]['counter']['p'] += 1
                else:
                    consensus[variantId]['counter']['v'] += 1
        return consensus

    def clear_tables(self):
        consensus = self.session.get('VKGL_consensus', num=10000)
        comments = self.session.get('VKGL_comments', num=10000)
        ids = []
        for row in consensus:
            ids.append(row['id'])
            self.old_diseases[row['id']] = ['' if 'disease' not in row else row['disease']['mim_number']][0]
        if len(ids) > 0:
            self.session.delete_list('VKGL_consensus', ids)
            print('Deleted consensus variants')
        ids = []
        for row in comments:
            id = row['id']
            if id.startswith('consensus_'):
                ids.append(id)
                self.old_comments[id] = row['comments']
        if len(ids) > 0:
            self.session.delete_list('VKGL_comments', ids)
            print('Deleted comments')

    def calculate_consensus(self):
        molgenis_table = []
        for id in self.lab_data:
            variant = self.lab_data[id]
            b = variant['counter']['b']
            p = variant['counter']['p']
            v = variant['counter']['v']
            if b > 1 and p == 0 and v == 0:
                variant['consensus_classification'] = '(Likely) benign (' + str(b) + ')'
            elif b == 0 and p > 1 and v == 0:
                variant['consensus_classification'] = '(Likely) pathogenic (' + str(p) + ')'
            elif b == 0 and p == 0 and v > 1:
                variant['consensus_classification'] = 'VUS(' + str(v) + ')'
            elif b > 0 and p > 0:
                variant['consensus_classification'] = 'Opposite classification'
            elif (b > 0 and v > 0) or (p > 0 and v > 0):
                variant['consensus_classification'] = 'No consensus'
            elif b == 1 or v == 1 or p == 1:
                variant['consensus_classification'] = 'Classified by one lab'
            self.lab_data[id] = variant
            del variant['counter']
            molgenis_table.append(variant)
        return molgenis_table

    def upload_comments(self):
        comments = []
        for id in self.lab_data:
            if id in self.old_comments:
                comments.append({'id': 'consensus_' + id, 'comments': self.old_comments[id]})
            else:
                comments.append({'id': 'consensus_' + id, 'comments': '-'})
        self.session.add_all('VKGL_comments', comments)

    def upload_consensus(self, entities):
        self.upload_comments()
        print('Comments uploaded')
        self.session.add_all('VKGL_consensus', entities)
        print('Consensus uploaded\nDone!')


def main():
    config = MolgenisConfigParser('config.txt').config
    labs = config['labs'].split(',')
    url = config['url']
    account = config['account']
    pwd = config['password']
    session = molgenis.Session(url)
    session.login(account, pwd)
    ConsensusTableGenerator(labs, session)


if __name__ == '__main__':
    main()
