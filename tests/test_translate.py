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
        expected = 'SXEWKTKRP'
        self.assertEqual(expected, seq_record.translate(table=1))

    def test_gapped_translation_with_mixed_codons(self):
        seq = 'TCTN--GAATGGAAGACAAAGCGTCCA'
        seq_record = SeqRecordExpanded(seq, reading_frame=1)
        result = seq_record.translate(table=1)
        self.assertEqual("SXEWKTKRP", result)

    def test_translation_with_missing_bp_as_question_mark(self):
        seq = "ATACGGTA?"
        seq_record = SeqRecordExpanded(seq, table=1, reading_frame=1,
                                       voucher_code="CP100-10", gene_code="wingless")
        result = seq_record.translate()
        self.assertEqual("IRX", result)

    def test_translation_with_missing_bp_as_dash(self):
        seq = 'ATACGGTA-'
        seq_record = SeqRecordExpanded(seq, table=1, reading_frame=1,
                                       voucher_code="CP100-10", gene_code="wingless")
        result = seq_record.translate()
        self.assertEqual("IRX", result)

    def test_translation_with_missing_bp_as_N(self):
        seq = 'ATACGGTAN'
        seq_record = SeqRecordExpanded(seq, table=1, reading_frame=1,
                                       voucher_code="CP100-10", gene_code="wingless")
        result = seq_record.translate()
        self.assertEqual("IRX", result)

    def test_translation_with_missing_bp_as_n(self):
        seq = 'ATACGGTAn'
        seq_record = SeqRecordExpanded(seq, table=1, reading_frame=1,
                                       voucher_code="CP100-10", gene_code="wingless")
        result = seq_record.translate()
        self.assertEqual("IRX", result)
