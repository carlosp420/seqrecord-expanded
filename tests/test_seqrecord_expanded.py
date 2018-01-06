import unittest

from degenerate_dna import exceptions

from seqrecord_expanded import SeqRecordExpanded


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
