import unittest

from seqrecord_expanded import SeqRecordExpanded


class TestCodonPositions(unittest.TestCase):
    def setUp(self):
        self.table = 1  # translation table
        self.voucher_code = 'CP100-09'
        self.taxonomy = {'genus': 'Aus', 'species': 'bus'}
        self.seq = 'GAATGGAAGACAAAGTCTCGTCCA'

    def test_missing_reading_frame(self):
        seq_record = SeqRecordExpanded(self.seq)
        self.assertRaises(AttributeError, seq_record.fist_codon_position)
        self.assertRaises(AttributeError, seq_record.second_codon_position)
        self.assertRaises(AttributeError, seq_record.third_codon_position)

    def test_wrong_reading_frame_int(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=4)
        self.assertRaises(ValueError, seq_record.fist_codon_position)
        self.assertRaises(ValueError, seq_record.second_codon_position)
        self.assertRaises(ValueError, seq_record.third_codon_position)

    def test_wrong_reading_frame_str(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame='1')
        self.assertRaises(ValueError, seq_record.fist_codon_position)
        self.assertRaises(ValueError, seq_record.second_codon_position)
        self.assertRaises(ValueError, seq_record.third_codon_position)

    def test_wrong_reading_frame_empty(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame='')
        self.assertRaises(AttributeError, seq_record.fist_codon_position)
        self.assertRaises(AttributeError, seq_record.second_codon_position)
        self.assertRaises(AttributeError, seq_record.third_codon_position)

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

    def test_getting_first_and_second_codon_positions_reading_frame_1(self):
        seq = 'GAATGGAAGACAAAGTCTCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = 'GATGAAACAATCCGCC'
        self.assertEqual(expected, seq_record.fist_and_second_codon_positions())
