from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from degenerate_dna import Degenera

from .utils import chain_and_flatten
from .exceptions import MissingParameterError


class SeqRecordExpanded(object):
    """Creates an Expanded SeqRecord.

    Assumes DNA ambiguous sequence.

    Arguments:
        - voucher_code       - code of voucher tha the sequence belongs to.
        - taxonomy           - dictionary {'genus': 'Aus', 'species': 'bus'}.
        - gene_code          - gene code.
        - reading_frame      - integer. 1, 2 or 3.
        - table  - integer. NCBI code for translation table.
    """
    def __init__(self, seq=None, voucher_code=None, taxonomy=None, gene_code=None,
                 reading_frame=None, table=None):
        self.seq = Seq(seq, alphabet=IUPAC.ambiguous_dna)
        self.voucher_code = voucher_code
        self.taxonomy = taxonomy
        self.gene_code = gene_code
        self.reading_frame = reading_frame
        self.table = table
        self._sequence_was_corrected = None

    def first_codon_position(self):
        """
        :return: string containing the first positions of each codon.
        """
        self._check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            first_position = seq[::3]
        elif self.reading_frame == 2:
            first_position = seq[1::3]
        else:  # self.reading_frame == 3
            first_position = seq[2::3]
        return first_position

    def _check_reading_frame(self):
        """
        Raises errors if reading frame is not integer and is not 1, 2 or 3.
        """
        if not self.reading_frame:
            raise AttributeError("The reading_frame attribute is not set.")
        if self.reading_frame not in [1, 2, 3]:
            raise ValueError("The reading_frame attribute should be either 1, 2 or 3.")

    def second_codon_position(self):
        """
        :return: string containing the second positions of each codon.
        """
        self._check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            second_position = seq[1::3]
        elif self.reading_frame == 2:
            second_position = seq[2::3]
        else:  # self.reading_frame == 3
            second_position = seq[::3]
        return second_position

    def third_codon_position(self):
        """
        :return: string containing the third positions of each codon.
        """
        self._check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            third_position = seq[2::3]
        elif self.reading_frame == 2:
            third_position = seq[::3]
        else:  # self.reading_frame == 3
            third_position = seq[1::3]
        return third_position

    def first_and_second_codon_positions(self):
        """
        :return: string containing both positions of each codon.
        """
        return chain_and_flatten(self.first_codon_position(), self.second_codon_position())

    def degenerate(self, method=None):
        self._check_reading_frame()
        self._correct_seq_based_on_reading_frame()

        if not method:
            table = self.table
            method = 'normal'
        else:
            table = 1
            method = method

        res = Degenera(dna=str(self.seq), table=table, method=method)
        res.degenerate()
        return res.degenerated

    def _correct_seq_based_on_reading_frame(self):
        """
        Trims leading end of `self.seq` if the reading frame does not start in 1st codon position
        of the sequence.
        """
        if not self._sequence_was_corrected:
            self._sequence_was_corrected = True

            if self.reading_frame == 2:
                self.seq = self.seq[1:]

            if self.reading_frame == 3:
                self.seq = self.seq[2:]

    def translate(self, table=None):
        """
        Uses BioPython translation method into Aminoacid sequence.
        :param table: Optional. It can be specified when creating the class instance.
        :return: str. Aminoacid sequence.
        """
        self._check_reading_frame()
        self._check_translation_table(table)
        self._correct_seq_based_on_reading_frame()

        if not table:
            return str(self.seq.translate(table=self.table))
        else:
            return str(self.seq.translate(table=table))

    def _check_translation_table(self, table):
        if self.table is None and table is None:
            raise MissingParameterError('It is necessary to specify the translation table to use:'
                                        ' seq_record.translate(table=1)')
