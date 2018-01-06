import unittest

from seqrecord_expanded import SeqRecordExpanded
from seqrecord_expanded.exceptions import MissingParameterError


class TestCodonPositions(unittest.TestCase):
    def setUp(self):
        self.table = 1  # translation table
        self.voucher_code = 'CP100-09'
        self.taxonomy = {'genus': 'A us-', 'species': 'bus(?)'}
        self.seq = 'GAATGGAAGACAAAGTCTCGTCCA'

    def test_fixing_taxon_name(self):
        seq_record = SeqRecordExpanded(
            self.seq,
            taxonomy=self.taxonomy,
            gene_code='wingless',
        )
        self.assertEqual("A_us_", seq_record.taxonomy['genus'])
        self.assertEqual("bus___", seq_record.taxonomy['species'])

    def test_missing_reading_frame(self):
        seq_record = SeqRecordExpanded(self.seq, gene_code='wingless')
        self.assertRaises(MissingParameterError, seq_record.first_codon_position)
        self.assertRaises(MissingParameterError, seq_record.second_codon_position)
        self.assertRaises(MissingParameterError, seq_record.third_codon_position)

    def test_wrong_reading_frame_int(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=4)
        self.assertRaises(ValueError, seq_record.first_codon_position)
        self.assertRaises(ValueError, seq_record.second_codon_position)
        self.assertRaises(ValueError, seq_record.third_codon_position)

    def test_wrong_reading_frame_str(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame='1')
        self.assertRaises(ValueError, seq_record.first_codon_position)
        self.assertRaises(ValueError, seq_record.second_codon_position)
        self.assertRaises(ValueError, seq_record.third_codon_position)

    def test_wrong_reading_frame_empty(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame='')
        self.assertRaises(ValueError, seq_record.first_codon_position)
        self.assertRaises(ValueError, seq_record.second_codon_position)
        self.assertRaises(ValueError, seq_record.third_codon_position)

    def test_getting_codon_positions_reading_frame_1(self):
        seq = '123123123123'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = '1111'
        self.assertEqual(expected, seq_record.first_codon_position(), 'Fist codon position')

        expected = '2222'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = '3333'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_codon_positions_reading_frame_2(self):
        seq = '23123123123'
        seq_record = SeqRecordExpanded(seq, reading_frame=2)
        expected = '111'
        self.assertEqual(expected, seq_record.first_codon_position(), 'Fist codon position')

        expected = '222'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = '333'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_codon_positions_reading_frame_3(self):
        seq = '3123123123123'
        seq_record = SeqRecordExpanded(seq, reading_frame=3)
        expected = '1111'
        self.assertEqual(expected, seq_record.first_codon_position(), 'Fist codon position')

        expected = '2222'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = '3333'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')


    def test_getting_first_and_second_codon_positions_reading_frame_1(self):
        seq = '123123123123'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = '12121212'
        self.assertEqual(expected, seq_record.first_and_second_codon_positions())

    def test_getting_first_and_second_codon_positions_reading_frame_2(self):
        seq = '23123123123123'
        seq_record = SeqRecordExpanded(seq, reading_frame=2)
        expected = '12121212'
        self.assertEqual(expected, seq_record.first_and_second_codon_positions())

    def test_getting_first_and_second_codon_positions_reading_frame_3(self):
        seq = '3123123123123'
        seq_record = SeqRecordExpanded(seq, reading_frame=3)
        expected = '12121212'
        self.assertEqual(expected, seq_record.first_and_second_codon_positions())
