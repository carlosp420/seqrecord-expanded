import unittest

from seqrecord_expanded import SeqRecordExpanded


class TestSeqRecordExpanded(unittest.TestCase):
    def setUp(self):
        self.seq = 'ACACGTCGACTCCGGCAAGTCCACTACCACAGGACATTTGATTTACAAATGTGGTGGTATCGACAAGCGT'
        self.gene_code = 'EF1a'
        self.reading_frame = 2
        self.table = 1  # translation table
        self.voucher_code = 'CP100-09'
        self.taxonomy = {'genus': 'Aus', 'species': 'bus'}

    def test_getting_fist_codon_position(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=self.reading_frame)
        result = seq_record.fist_codon_position()
        expected = 'CGGTGATAAAGCTATATGGAGAC'
        self.assertEqual(expected, result)

    def test_getting_second_codon_position(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=self.reading_frame)
        result = seq_record.second_codon_position()
        expected = 'ATACGACCCCGATTAAGGGTAAG'
        self.assertEqual(expected, result)
