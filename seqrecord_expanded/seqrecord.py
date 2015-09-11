from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from .utils import chain_and_flatten


class SeqRecordExpanded(SeqRecord):
    """Creates an Expanded SeqRecord.

    Assumes DNA ambiguous sequence.

    Arguments:
        - voucher_code       - code of voucher tha the sequence belongs to.
        - taxonomy           - dictionary {'genus': 'Aus', 'species': 'bus'}.
        - gene_code          - gene code.
        - reading_frame      - integer. 1, 2 or 3.
        - translation_table  - integer. NCBI code for translation table.
    """
    def __init__(self, *args, voucher_code=None, taxonomy=None, gene_code=None,
                 reading_frame=None, translation_table=None, **kwargs):
        super(SeqRecordExpanded, self).__init__(*args, **kwargs)
        self._seq = Seq(args[0], alphabet=IUPAC.ambiguous_dna)
        self.voucher_code = voucher_code
        self.taxonomy = taxonomy
        self.gene_code = gene_code
        self.reading_frame = reading_frame
        self.translation_table = translation_table

    def fist_codon_position(self):
        """
        :return: string containing the first positions of each codon.
        """
        self.check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            first_position = seq[::3]
        elif self.reading_frame == 2:
            first_position = seq[1::3]
        else:  # self.reading_frame == 3
            first_position = seq[2::3]
        return first_position

    def check_reading_frame(self):
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
        self.check_reading_frame()

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
        self.check_reading_frame()

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
        return chain_and_flatten(self.fist_codon_position(), self.second_codon_position())
