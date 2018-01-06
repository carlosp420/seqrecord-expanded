import re
import warnings

from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from degenerate_dna import Degenera

from .utils import chain_and_flatten
from .exceptions import MissingParameterError, TranslationErrorMixedGappedSeq
from ._warnings import SeqRecordExpandedWarning


class SeqRecordExpanded(object):
    """Creates an Expanded SeqRecord.

    Assumes DNA ambiguous sequence.

    Parameters:
        seq (str):            DNA sequence
        voucher_code (str):   code of voucher that the sequence belongs to
        taxonomy (dict):      ``{'genus': 'Aus', 'species': 'bus'}``
        lineage (str):        ``Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta;
                              Pterygota; Neoptera; Holometabola; Lepidoptera; Glossata; Ditrysia;
                              Papilionoidea; Nymphalidae; Satyrinae; Satyrini; Euptychiina;``
        gene_code (str):      gene code
        reading_frame (int):  1, 2 or 3.
        table (int):          NCBI code for translation table

    Attributes:
        seq:               DNA sequence as string
        voucher_code:      Code of voucher tha the sequence belongs to.
        taxonomy:          Dictionary ``{'genus': 'Aus', 'species': 'bus'}``.
        gene_code:         Gene code.
        reading_frame:     1, 2 or 3.
        table:             NCBI code for translation table.
        warnings:          List.

    Raises:
        MissingParameterError:  if user wants either first, second, third or
                                first and second codon positions and
                                ``reading_frame`` is not specified.

    """
    def __init__(self, seq=None, voucher_code=None, taxonomy=None, lineage=None,
                 gene_code=None, reading_frame=None, table=None):
        self.warnings = []
        self.seq = Seq(
            seq.replace("-", "?"),
            alphabet=IUPAC.ambiguous_dna,
        )
        self.voucher_code = voucher_code
        self.taxonomy = ""
        self.lineage = lineage
        self.gene_code = gene_code
        self.reading_frame = reading_frame
        self.table = table
        self._sequence_was_corrected = None
        self._clean_taxonomy(taxonomy)

    def _clean_taxonomy(self, taxonomy):
        self.taxonomy = dict()
        if taxonomy:
            for key, value in taxonomy.items():
                # remove special characters so Biopython will not choke on them.
                self.taxonomy[key] = re.sub("\W", "_", value)

    def first_codon_position(self):
        """
        :return: string containing the first positions of each codon.

        """
        self._check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            pass
        elif self.reading_frame == 2:
            seq = seq[2:]
        elif self.reading_frame == 3:
            seq = seq[1:]
        else:  # None
            raise MissingParameterError('reading_frame attribute for gene {0} '
                                        'should be either 1, 2 or 3.'.format(self.gene_code))
        return seq[::3]

    def _check_reading_frame(self):
        """Raises errors if reading frame is not integer and is not 1, 2, 3 or None.

        """
        if self.reading_frame not in [1, 2, 3, None]:
            raise ValueError("The reading_frame attribute should be either 1, 2, 3 or None.")

    def second_codon_position(self):
        """
        :return: string containing the second positions of each codon.

        """
        self._check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            pass
        elif self.reading_frame == 2:
            seq = seq[2:]
        elif self.reading_frame == 3:
            seq = seq[1:]
        else:  # None
            raise MissingParameterError('reading_frame attribute for gene {0} '
                                        'should be either 1, 2 or 3.'.format(self.gene_code))
        return seq[1::3]

    def third_codon_position(self):
        """
        :return: string containing the third positions of each codon.

        """
        self._check_reading_frame()

        seq = str(self.seq)
        if self.reading_frame == 1:
            pass
        elif self.reading_frame == 2:
            seq = seq[2:]
        elif self.reading_frame == 3:
            seq = seq[1:]
        else:  # None
            raise MissingParameterError('reading_frame attribute for gene {0} '
                                        'should be either 1, 2 or 3.'.format(self.gene_code))
        return seq[2::3]

    def first_and_second_codon_positions(self):
        """
        :return: string containing both positions of each codon.

        """
        return chain_and_flatten(self.first_codon_position(), self.second_codon_position())

    def degenerate(self, method=None):
        """
        Parameters:
            method (str):   S, Z, SZ, normal

        Returns:
            (str): Degenerated sequence using Zwick et al methods.

        """
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
        """Trims leading end of `self.seq`.

        If the reading frame does not start in 1st codon position of the
        sequence.

        """
        if not self._sequence_was_corrected:
            self._sequence_was_corrected = True

            if self.reading_frame == 1:
                pass
            elif self.reading_frame == 2:
                self.seq = self.seq[1:]
            elif self.reading_frame == 3:
                self.seq = self.seq[2:]
            else:  # reading_frame is None
                self.seq = '?'
                msg = 'reading_frame attribute should be either 1, 2 or 3.'
                warnings.warn(msg, SeqRecordExpandedWarning)
                self.warnings.append(msg)

    def translate(self, table=None):
        """Uses BioPython translation method into Aminoacid sequence.

        Parameters:
            table (int): Optional. It can be specified when creating the class instance.

        Returns:
            (str): Aminoacid sequence.

        """
        if not table:
            table = self.table

        self._check_reading_frame()
        self._check_translation_table(table)
        self._correct_seq_based_on_reading_frame()

        new_seq = Seq(str(self.seq).replace('?', 'N'), alphabet=IUPAC.ambiguous_dna)
        try:
            translated_seq = self._translate(new_seq, table)
        except TranslationError as e:
            raise TranslationErrorMixedGappedSeq(self.voucher_code, self.gene_code, e)
        return translated_seq

    def _translate(self, new_seq, table):
        if not table:
            return str(new_seq.translate(table=self.table, gap="-"))
        else:
            return str(new_seq.translate(table=table, gap="-"))

    def _check_translation_table(self, table):
        if self.table is None and table is None:
            raise MissingParameterError('It is necessary to specify the translation'
                                        ' table to use: seq_record.translate(table=1)')
