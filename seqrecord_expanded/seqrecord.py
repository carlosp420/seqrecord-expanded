from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


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
