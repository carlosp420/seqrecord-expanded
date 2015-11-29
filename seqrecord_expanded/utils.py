import itertools
import warnings

import six
if six.PY2:
    from six.moves import zip_longest
else:
    from itertools import zip_longest

from Bio._py3k import range
from Bio._py3k import basestring
from Bio import Alphabet
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


def chain_and_flatten(seq1, seq2):
    """Takes two strings (first and second codon positions) and chains and
    intercalates them.

    Returns:
        (str): String of combining the two seq strings.
    """
    my_chain = zip_longest(seq1, seq2)
    out = [i for i in itertools.chain.from_iterable(my_chain) if i]
    return ''.join(out)


# Implementation of gapped translation. See pull request in Biopython
# https://github.com/biopython/biopython/pull/661
# Remove once Biopython implement this.
class NewSeq(Seq):
    """A read-only sequence object (essentially a string with an alphabet).

    Like normal python strings, our basic sequence object is immutable.
    This prevents you from doing my_seq[5] = "A" for example, but does allow
    Seq objects to be used as dictionary keys.

    The Seq object provides a number of string like methods (such as count,
    find, split and strip), which are alphabet aware where appropriate.

    In addition to the string like sequence, the Seq object has an alphabet
    property. This is an instance of an Alphabet class from Bio.Alphabet,
    for example generic DNA, or IUPAC DNA. This describes the type of molecule
    (e.g. RNA, DNA, protein) and may also indicate the expected symbols
    (letters).

    The Seq object also provides some biological methods, such as complement,
    reverse_complement, transcribe, back_transcribe and translate (which are
    not applicable to sequences with a protein alphabet).
    """
    def translate(self, table="Standard", stop_symbol="*", to_stop=False,
                  cds=False, gap=None):
        """Turns a nucleotide sequence into a protein sequence. New Seq object.

        This method will translate DNA or RNA sequences, and those with a
        nucleotide or generic alphabet.  Trying to translate a protein
        sequence raises an exception.

        Arguments:
            - table - Which codon table to use?  This can be either a name
              (string), an NCBI identifier (integer), or a CodonTable
              object (useful for non-standard genetic codes).  This
              defaults to the "Standard" table.
            - stop_symbol - Single character string, what to use for terminators.
              This defaults to the asterisk, "*".
            - to_stop - Boolean, defaults to False meaning do a full translation
              continuing on past any stop codons (translated as the
              specified stop_symbol).  If True, translation is
              terminated at the first in frame stop codon (and the
              stop_symbol is not appended to the returned protein
              sequence).
            - cds - Boolean, indicates this is a complete CDS.  If True,
              this checks the sequence starts with a valid alternative start
              codon (which will be translated as methionine, M), that the
              sequence length is a multiple of three, and that there is a
              single in frame stop codon at the end (this will be excluded
              from the protein sequence, regardless of the to_stop option).
              If these tests fail, an exception is raised.
            - gap - Single character string to denote symbol used for gaps.
              It will try to guess the gap character from the alphabet.

        e.g. Using the standard table:

        >>> coding_dna = Seq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> coding_dna.translate()
        Seq('VAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> coding_dna.translate(stop_symbol="@")
        Seq('VAIVMGR@KGAR@', HasStopCodon(ExtendedIUPACProtein(), '@'))
        >>> coding_dna.translate(to_stop=True)
        Seq('VAIVMGR', ExtendedIUPACProtein())

        Now using NCBI table 2, where TGA is not a stop codon:

        >>> coding_dna.translate(table=2)
        Seq('VAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> coding_dna.translate(table=2, to_stop=True)
        Seq('VAIVMGRWKGAR', ExtendedIUPACProtein())

        In fact, GTG is an alternative start codon under NCBI table 2, meaning
        this sequence could be a complete CDS:

        >>> coding_dna.translate(table=2, cds=True)
        Seq('MAIVMGRWKGAR', ExtendedIUPACProtein())

        It isn't a valid CDS under NCBI table 1, due to both the start codon and
        also the in frame stop codons:

        >>> coding_dna.translate(table=1, cds=True)
        Traceback (most recent call last):
            ...
        TranslationError: First codon 'GTG' is not a start codon

        If the sequence has no in-frame stop codon, then the to_stop argument
        has no effect:

        >>> coding_dna2 = Seq("TTGGCCATTGTAATGGGCCGC")
        >>> coding_dna2.translate()
        Seq('LAIVMGR', ExtendedIUPACProtein())
        >>> coding_dna2.translate(to_stop=True)
        Seq('LAIVMGR', ExtendedIUPACProtein())

        When translating gapped sequences, the gap character is inferred from
        the alphabet:

        >>> from Bio.Alphabet import Gapped
        >>> coding_dna3 = Seq("GTG---GCCATT", Gapped(IUPAC.unambiguous_dna))
        >>> coding_dna3.translate()
        Seq('V-AI', Gapped(ExtendedIUPACProtein(), '-'))

        It is possible to pass the gap character when the alphabet is missing:

        >>> coding_dna4 = Seq("GTG---GCCATT")
        >>> coding_dna4.translate(gap='-')
        Seq('V-AI', Gapped(ExtendedIUPACProtein(), '-'))

        NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
        or a stop codon.  These are translated as "X".  Any invalid codon
        (e.g. "TA?" or "T-A") will throw a TranslationError.

        NOTE - This does NOT behave like the python string's translate
        method.  For that use str(my_seq).translate(...) instead.
        """
        if isinstance(table, str) and len(table) == 256:
            raise ValueError("The Seq object translate method DOES NOT take "
                             "a 256 character string mapping table like "
                             "the python string object's translate method. "
                             "Use str(my_seq).translate(...) instead.")
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet):
            raise ValueError("Proteins cannot be translated!")
        try:
            table_id = int(table)
        except ValueError:
            # Assume its a table name
            if self.alphabet == IUPAC.unambiguous_dna:
                # Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_dna_by_name[table]
            elif self.alphabet == IUPAC.unambiguous_rna:
                # Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_rna_by_name[table]
            else:
                # This will use the extended IUPAC protein alphabet with X etc.
                # The same table can be used for RNA or DNA (we use this for
                # translating strings).
                codon_table = CodonTable.ambiguous_generic_by_name[table]
        except (AttributeError, TypeError):
            # Assume its a CodonTable object
            if isinstance(table, CodonTable.CodonTable):
                codon_table = table
            else:
                raise ValueError('Bad table argument')
        else:
            # Assume its a table ID
            if self.alphabet == IUPAC.unambiguous_dna:
                # Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_dna_by_id[table_id]
            elif self.alphabet == IUPAC.unambiguous_rna:
                # Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_rna_by_id[table_id]
            else:
                # This will use the extended IUPAC protein alphabet with X etc.
                # The same table can be used for RNA or DNA (we use this for
                # translating strings).
                codon_table = CodonTable.ambiguous_generic_by_id[table_id]

        # Deal with gaps for translation
        if hasattr(self.alphabet, "gap_char"):
            if not gap:
                gap = self.alphabet.gap_char
            elif gap != self.alphabet.gap_char:
                raise ValueError("Gap {0!r} does not match {1!r} from alphabet".format(
                    gap, self.alphabet.gap_char))

        protein = _translate_str(str(self), codon_table, stop_symbol, to_stop,
                                 cds, gap=gap)

        if gap and gap in protein:
            alphabet = Alphabet.Gapped(codon_table.protein_alphabet, gap)
        else:
            alphabet = codon_table.protein_alphabet

        if stop_symbol in protein:
            alphabet = Alphabet.HasStopCodon(alphabet, stop_symbol)

        return Seq(protein, alphabet)


def _translate_str(sequence, table, stop_symbol="*", to_stop=False,
                   cds=False, pos_stop="X", gap=None):
    """Helper function to translate a nucleotide string (PRIVATE).

    Arguments:
        - sequence - a string
        - table - a CodonTable object (NOT a table name or id number)
        - stop_symbol - a single character string, what to use for terminators.
        - to_stop - boolean, should translation terminate at the first
          in frame stop codon?  If there is no in-frame stop codon
          then translation continues to the end.
        - pos_stop - a single character string for a possible stop codon
          (e.g. TAN or NNN)
        - cds - Boolean, indicates this is a complete CDS.  If True, this
          checks the sequence starts with a valid alternative start
          codon (which will be translated as methionine, M), that the
          sequence length is a multiple of three, and that there is a
          single in frame stop codon at the end (this will be excluded
          from the protein sequence, regardless of the to_stop option).
          If these tests fail, an exception is raised.
        - gap - Single character string to denote symbol used for gaps.
          Defaults to None.

    Returns a string.

    e.g.

    >>> from Bio.Data import CodonTable
    >>> table = CodonTable.ambiguous_dna_by_id[1]
    >>> _translate_str("AAA", table)
    'K'
    >>> _translate_str("TAR", table)
    '*'
    >>> _translate_str("TAN", table)
    'X'
    >>> _translate_str("TAN", table, pos_stop="@")
    '@'
    >>> _translate_str("TA?", table)
    Traceback (most recent call last):
       ...
    TranslationError: Codon 'TA?' is invalid

    In a change to older versions of Biopython, partial codons are now
    always regarded as an error (previously only checked if cds=True)
    and will trigger a warning (likely to become an exception in a
    future release).

    If **cds=True**, the start and stop codons are checked, and the start
    codon will be translated at methionine. The sequence must be an
    while number of codons.

    >>> _translate_str("ATGCCCTAG", table, cds=True)
    'MP'
    >>> _translate_str("AAACCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    TranslationError: First codon 'AAA' is not a start codon
    >>> _translate_str("ATGCCCTAGCCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    TranslationError: Extra in frame stop codon found.
    """
    sequence = sequence.upper()
    amino_acids = []
    forward_table = table.forward_table
    stop_codons = table.stop_codons
    if table.nucleotide_alphabet.letters is not None:
        valid_letters = set(table.nucleotide_alphabet.letters.upper())
    else:
        # Assume the worst case, ambiguous DNA or RNA:
        valid_letters = set(IUPAC.ambiguous_dna.letters.upper() +
                            IUPAC.ambiguous_rna.letters.upper())
    n = len(sequence)
    if cds:
        if str(sequence[:3]).upper() not in table.start_codons:
            raise CodonTable.TranslationError(
                "First codon '{0}' is not a start codon".format(sequence[:3]))
        if n % 3 != 0:
            raise CodonTable.TranslationError(
                "Sequence length {0} is not a multiple of three".format(n))
        if str(sequence[-3:]).upper() not in stop_codons:
            raise CodonTable.TranslationError(
                "Final codon '{0}' is not a stop codon".format(sequence[-3:]))
        # Don't translate the stop symbol, and manually translate the M
        sequence = sequence[3:-3]
        n -= 6
        amino_acids = ["M"]
    elif n % 3 != 0:
        from Bio import BiopythonWarning
        warnings.warn("Partial codon, len(sequence) not a multiple of three. "
                      "Explicitly trim the sequence or add trailing N before "
                      "translation. This may become an error in future.",
                      BiopythonWarning)
    if gap is not None:
        if not isinstance(gap, basestring):
            raise TypeError("Gap character should be a single character string.")
        elif len(gap) > 1:
            raise ValueError("Gap character should be a single character string.")

    for i in range(0, n - n % 3, 3):
        codon = sequence[i:i + 3]
        try:
            amino_acids.append(forward_table[codon])
        except (KeyError, CodonTable.TranslationError):
            if codon in table.stop_codons:
                if cds:
                    raise CodonTable.TranslationError(
                        "Extra in frame stop codon found.")
                if to_stop:
                    break
                amino_acids.append(stop_symbol)
            elif valid_letters.issuperset(set(codon)):
                # Possible stop codon (e.g. NNN or TAN)
                amino_acids.append(pos_stop)
            elif gap is not None and codon == gap * 3:
                # Gapped translation
                amino_acids.append(gap)
            else:
                raise CodonTable.TranslationError(
                    "Codon '{0}' is invalid".format(codon))
    return "".join(amino_acids)


def translate(sequence, table="Standard", stop_symbol="*", to_stop=False,
              cds=False, gap=None):
    """Translate a nucleotide sequence into amino acids.

    If given a string, returns a new string object. Given a Seq or
    MutableSeq, returns a Seq object with a protein alphabet.

    Arguments:
        - table - Which codon table to use?  This can be either a name (string),
          an NCBI identifier (integer), or a CodonTable object (useful
          for non-standard genetic codes).  Defaults to the "Standard"
          table.
        - stop_symbol - Single character string, what to use for any
          terminators, defaults to the asterisk, "*".
        - to_stop - Boolean, defaults to False meaning do a full
          translation continuing on past any stop codons
          (translated as the specified stop_symbol).  If
          True, translation is terminated at the first in
          frame stop codon (and the stop_symbol is not
          appended to the returned protein sequence).
        - cds - Boolean, indicates this is a complete CDS.  If True, this
          checks the sequence starts with a valid alternative start
          codon (which will be translated as methionine, M), that the
          sequence length is a multiple of three, and that there is a
          single in frame stop codon at the end (this will be excluded
          from the protein sequence, regardless of the to_stop option).
          If these tests fail, an exception is raised.
        - gap - Single character string to denote symbol used for gaps.
          Defaults to None.

    A simple string example using the default (standard) genetic code:

    >>> coding_dna = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    >>> translate(coding_dna)
    'VAIVMGR*KGAR*'
    >>> translate(coding_dna, stop_symbol="@")
    'VAIVMGR@KGAR@'
    >>> translate(coding_dna, to_stop=True)
    'VAIVMGR'

    Now using NCBI table 2, where TGA is not a stop codon:

    >>> translate(coding_dna, table=2)
    'VAIVMGRWKGAR*'
    >>> translate(coding_dna, table=2, to_stop=True)
    'VAIVMGRWKGAR'

    In fact this example uses an alternative start codon valid under NCBI table 2,
    GTG, which means this example is a complete valid CDS which when translated
    should really start with methionine (not valine):

    >>> translate(coding_dna, table=2, cds=True)
    'MAIVMGRWKGAR'

    Note that if the sequence has no in-frame stop codon, then the to_stop
    argument has no effect:

    >>> coding_dna2 = "GTGGCCATTGTAATGGGCCGC"
    >>> translate(coding_dna2)
    'VAIVMGR'
    >>> translate(coding_dna2, to_stop=True)
    'VAIVMGR'

    NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
    or a stop codon.  These are translated as "X".  Any invalid codon
    (e.g. "TA?" or "T-A") will throw a TranslationError.

    It will however translate either DNA or RNA.
    """
    if isinstance(sequence, Seq):
        return sequence.translate(table, stop_symbol, to_stop, cds)
    elif isinstance(sequence, MutableSeq):
        # Return a Seq object
        return sequence.toseq().translate(table, stop_symbol, to_stop, cds)
    else:
        # Assume its a string, return a string
        try:
            codon_table = CodonTable.ambiguous_generic_by_id[int(table)]
        except ValueError:
            codon_table = CodonTable.ambiguous_generic_by_name[table]
        except (AttributeError, TypeError):
            if isinstance(table, CodonTable.CodonTable):
                codon_table = table
            else:
                raise ValueError('Bad table argument')
        return _translate_str(sequence, codon_table, stop_symbol, to_stop, cds,
                              gap=gap)
