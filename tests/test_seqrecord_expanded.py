import unittest

from degenerate_dna import exceptions

from seqrecord_expanded import SeqRecordExpanded
from seqrecord_expanded.exceptions import MissingParameterError


class TestCodonPositions(unittest.TestCase):
    def setUp(self):
        self.table = 1  # translation table
        self.voucher_code = 'CP100-09'
        self.taxonomy = {'genus': 'Aus', 'species': 'bus'}
        self.seq = 'GAATGGAAGACAAAGTCTCGTCCA'

    def test_missing_reading_frame(self):
        seq_record = SeqRecordExpanded(self.seq)
        self.assertEqual('?', seq_record.first_codon_position())
        self.assertEqual(seq_record.warnings, ['reading_frame attribute should be either 1, 2 or 3.'])
        self.assertEqual('?', seq_record.second_codon_position())
        self.assertEqual('?', seq_record.third_codon_position())

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
        seq = 'GAATGGAAGACAAAGTCTCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = 'GTAAATCC'
        self.assertEqual(expected, seq_record.first_codon_position(), 'Fist codon position')

        expected = 'AGACACGC'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = 'AGGAGTTA'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_codon_positions_reading_frame_2(self):
        seq = 'ACACGTCGACTCCGGCAAGTCCACTACCACAGGACATTTGATTTACAAATGTGGTGGTATCGACAAGCGT'
        seq_record = SeqRecordExpanded(seq, reading_frame=2)
        expected = 'CGGTGATAAAGCTATATGGAGAC'
        self.assertEqual(expected, seq_record.first_codon_position(), 'Fist codon position')

        expected = 'ATACGACCCCGATTAAGGGTAAG'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = 'ACCCCCGCTCAATGTCATTTCCGT'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_codon_positions_reading_frame_3(self):
        seq = 'GTCGTGGGGGCCCACGTGGACGTGG'
        seq_record = SeqRecordExpanded(seq, reading_frame=3)
        expected = 'CGGCCGCG'
        self.assertEqual(expected, seq_record.first_codon_position(), 'Fist codon position')

        expected = 'GGGGCGGGG'
        self.assertEqual(expected, seq_record.second_codon_position(), 'Second codon position')

        expected = 'TTGCATAT'
        self.assertEqual(expected, seq_record.third_codon_position(), 'Third codon position')

    def test_getting_first_and_second_codon_positions_reading_frame_1(self):
        seq = 'GAATGGAAGACAAAGTCTCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = 'GATGAAACAATCCGCC'
        self.assertEqual(expected, seq_record.first_and_second_codon_positions())


class TestDegenerate(unittest.TestCase):
    def setUp(self):
        self.table = 1  # translation table
        self.voucher_code = 'CP100-09'
        self.taxonomy = {'genus': 'Aus', 'species': 'bus'}
        self.seq = 'TCTGAATGGAAGACAAAGCGTCCA'

    def test_degen_no_reading_frame(self):
        seq_record = SeqRecordExpanded(self.seq)
        self.assertRaises(exceptions.MissingParameterError, seq_record.degenerate)

    def test_degen_missing_table_and_method(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=1)
        self.assertRaises(exceptions.WrongParameterError, seq_record.degenerate, 'Missing reading_frame.')

    def test_degen_standard(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=1, table=1)
        expected = 'TCNGARTGGAARACNAARMGNCCN'
        self.assertEqual(expected, seq_record.degenerate(), 'Using reading_frame=1')

        seq_record = SeqRecordExpanded(self.seq, reading_frame=2, table=1)
        expected = 'YTNAAYGGNMGNCARAGYGTNCA'
        self.assertEqual(expected, seq_record.degenerate(), 'Using reading_frame=2')

    def test_degen_s(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=1)
        expected = 'AGYGARTGGAARACNAARMGNCCN'
        self.assertEqual(expected, seq_record.degenerate(method='S'), 'Using reading_frame=1')

        seq_record = SeqRecordExpanded(self.seq, reading_frame=3)
        expected = 'TGAATGGARGAYAARGCNAGYA'
        self.assertEqual(expected, seq_record.degenerate(method='S'), 'Using reading_frame=3')

    def test_degen_z(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=1)
        expected = 'TCNGARTGGAARACNAARMGNCCN'
        self.assertEqual(expected, seq_record.degenerate(method='Z'))

    def test_degen_sz(self):
        seq_record = SeqRecordExpanded(self.seq, reading_frame=1)
        expected = 'NNNGARTGGAARACNAARMGNCCN'
        self.assertEqual(expected, seq_record.degenerate(method='SZ'), 'Using reading_frame=1')

        seq_record = SeqRecordExpanded(self.seq, reading_frame=3)
        expected = 'TGAATGGARGAYAARGCNNNNA'
        self.assertEqual(expected, seq_record.degenerate(method='SZ'), 'Using reading_frame=3')


class TestTranslate(unittest.TestCase):
    def setUp(self):
        pass

    def test_translate(self):
        seq = 'TCTGAATGGAAGACAAAGCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1, table=1)
        expected = 'SEWKTKRP'
        self.assertEqual(expected, seq_record.translate(), 'Using reading_frame=1')

        seq = 'TCTGAATGGAA?ACAAAGCGT???'
        seq_record = SeqRecordExpanded(seq, reading_frame=1, table=1)
        expected = 'SEWXTKRX'
        self.assertEqual(expected, seq_record.translate(), 'Using reading_frame=1')

        seq = 'ACACGTCGACTCCGGCAAGTCCACTACCACAGGA'
        seq_record = SeqRecordExpanded(seq, reading_frame=2, table=1)
        expected = 'HVDSGKSTTTG'
        self.assertEqual(expected, seq_record.translate(), 'Using reading_frame=2')

    def test_translate_no_table_at_creation_class_instance(self):
        seq = 'TCTGAATGGAAGACAAAGCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        self.assertRaises(MissingParameterError, seq_record.translate)

        seq = 'ACACGTCGACTCCGGCAAGTCCACTACCACAGGA'
        seq_record = SeqRecordExpanded(seq, reading_frame=2)
        self.assertRaises(MissingParameterError, seq_record.translate)

    def test_translate_with_table_at_function_level(self):
        seq = 'TCTGAATGGAAGACAAAGCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = 'SEWKTKRP'
        self.assertEqual(expected, seq_record.translate(table=1))

        seq = 'ACACGTCGACTCCGGCAAGTCCACTACCACAGGA'
        seq_record = SeqRecordExpanded(seq, reading_frame=2)
        expected = 'HVDSGKSTTTG'
        self.assertEqual(expected, seq_record.translate(table=1))
