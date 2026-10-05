from collections import namedtuple
import unittest
from pathlib import Path
from spliceai.utils import Annotator, get_delta_scores, INFO_FIELD_KEYS


Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])


def info_fields(scores):
    # the SpliceAI INFO entries the CLI writes to the VCF for these scores
    return ['|'.join(str(d[key]) for key in INFO_FIELD_KEYS) for d in scores]


class TestDeltaScore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        data_dir = Path(__file__).parent / 'data'
        fasta_path = str(data_dir / 'test.fa')
        fasta_without_prefix_path = str(data_dir / 'test_without_prefix.fa')
        cls.ann = Annotator(fasta_path, 'grch37')
        cls.ann_without_prefix = Annotator(fasta_without_prefix_path, 'grch37')

    def test_get_delta_score_acceptor(self):

        record = Record('10', 94077, 'A', ['C'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(info_fields(scores), ['C|TUBB8|-|0.147|0.267|0.000|0.051|89|-23|-267|193|0.524|0.671|0.979|0.712|0.001|0.001|0.133|0.082'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(info_fields(scores), ['C|TUBB8|-|0.147|0.267|0.000|0.051|89|-23|-267|193|0.524|0.671|0.979|0.712|0.001|0.001|0.133|0.082'])

        record = Record('chr10', 94077, 'A', ['C'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(info_fields(scores), ['C|TUBB8|-|0.147|0.267|0.000|0.051|89|-23|-267|193|0.524|0.671|0.979|0.712|0.001|0.001|0.133|0.082'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(info_fields(scores), ['C|TUBB8|-|0.147|0.267|0.000|0.051|89|-23|-267|193|0.524|0.671|0.979|0.712|0.001|0.001|0.133|0.082'])

    def test_get_delta_score_donor(self):

        record = Record('10', 94555, 'C', ['T'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(info_fields(scores), ['T|TUBB8|-|0.006|0.182|0.154|0.622|-2|110|-190|0|0.000|0.006|0.974|0.792|0.556|0.710|0.986|0.364'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(info_fields(scores), ['T|TUBB8|-|0.006|0.182|0.154|0.622|-2|110|-190|0|0.000|0.006|0.974|0.792|0.556|0.710|0.986|0.364'])

        record = Record('chr10', 94555, 'C', ['T'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(info_fields(scores), ['T|TUBB8|-|0.006|0.182|0.154|0.622|-2|110|-190|0|0.000|0.006|0.974|0.792|0.556|0.710|0.986|0.364'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(info_fields(scores), ['T|TUBB8|-|0.006|0.182|0.154|0.622|-2|110|-190|0|0.000|0.006|0.974|0.792|0.556|0.710|0.986|0.364'])

    def test_get_delta_score_mnp(self):

        record = Record('10', 94077, 'ACT', ['CCT'])
        scores = get_delta_scores(record, self.ann, 50, 0)
        self.assertEqual(info_fields(scores), ['CCT|TUBB8|-|0.073|0.267|0.000|0.015|0|-23|19|-22|0.027|0.101|0.979|0.712|0.000|0.000|0.105|0.090'])
        scores = get_delta_scores(record, self.ann_without_prefix, 50, 0)
        self.assertEqual(info_fields(scores), ['CCT|TUBB8|-|0.073|0.267|0.000|0.015|0|-23|19|-22|0.027|0.101|0.979|0.712|0.000|0.000|0.105|0.090'])

        record = Record('10', 94555, 'CGA', ['TGA'])
        scores = get_delta_scores(record, self.ann, 50, 0)
        self.assertEqual(info_fields(scores), ['TGA|TUBB8|-|0.006|0.000|0.112|0.622|-2|-6|-23|0|0.000|0.006|0.000|0.000|0.006|0.118|0.986|0.364'])
        scores = get_delta_scores(record, self.ann_without_prefix, 50, 0)
        self.assertEqual(info_fields(scores), ['TGA|TUBB8|-|0.006|0.000|0.112|0.622|-2|-6|-23|0|0.000|0.006|0.000|0.000|0.006|0.118|0.986|0.364'])
