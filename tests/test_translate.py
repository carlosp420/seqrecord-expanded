import unittest

from degenerate_dna import exceptions

from seqrecord_expanded import SeqRecordExpanded
from seqrecord_expanded.exceptions import MissingParameterError
from seqrecord_expanded.exceptions import TranslationErrorMixedGappedSeq


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

    def test_gapped_translation(self):
        seq = 'TCT---GAATGGAAGACAAAGCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        expected = 'S-EWKTKRP'
        self.assertEqual(expected, seq_record.translate(table=1))

    def test_gapped_translation_with_mixed_codons(self):
        seq = 'TCTN--GAATGGAAGACAAAGCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        self.assertRaises(TranslationErrorMixedGappedSeq, seq_record.translate,
                          table=1)

        try:
            seq_record.translate(table=1)
        except TranslationErrorMixedGappedSeq as e:
            self.assertTrue("Gene" in e.__str__())
