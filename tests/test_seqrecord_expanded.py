import unittest

from seqrecord_expanded import SeqRecordExpanded


class TestCodonPositions(unittest.TestCase):
    def setUp(self):
        self.table = 1  # translation table
        self.voucher_code = 'CP100-09'
        self.taxonomy = {'genus': 'Aus', 'species': 'bus'}

    def test_getting_codon_positions_reading_frame_1(self):
        seq = 'GAATGGAAGACAAAGTCTCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = 'GTAAATCC'
        self.assertEqual(expected, seq_record.fist_codon_position(), 'Fist codon position')

        expected = 'AGACACGC'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = 'AGGAGTTA'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_codon_positions_reading_frame_2(self):
        seq = 'ACACGTCGACTCCGGCAAGTCCACTACCACAGGACATTTGATTTACAAATGTGGTGGTATCGACAAGCGT'
        seq_record = SeqRecordExpanded(seq, reading_frame=2)
        expected = 'CGGTGATAAAGCTATATGGAGAC'
        self.assertEqual(expected, seq_record.fist_codon_position(), 'Fist codon position')

        expected = 'ATACGACCCCGATTAAGGGTAAG'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = 'ACCCCCGCTCAATGTCATTTCCGT'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_codon_positions_reading_frame_3(self):
        seq = 'GTCGTGGGGGCCCACGTGGACGTGG'
        seq_record = SeqRecordExpanded(seq, reading_frame=3)
        expected = 'CGGCCGCG'
        self.assertEqual(expected, seq_record.fist_codon_position(), 'Fist codon position')

        expected = 'GGGGCGGGG'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = 'TTGCATAT'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')
