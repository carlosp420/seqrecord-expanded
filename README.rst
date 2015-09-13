.. image:: https://img.shields.io/pypi/v/seqrecord-expanded.svg
   :target: https://pypi.python.org/pypi/seqrecord-expanded/
   :alt: Latest PyPI version

.. image:: https://travis-ci.org/carlosp420/seqrecord-expanded.png
   :target: https://travis-ci.org/carlosp420/seqrecord-expanded
   :alt: Latest Travis CI build status

.. image:: https://coveralls.io/repos/carlosp420/seqrecord-expanded/badge.svg?branch=master&service=github
   :target: https://coveralls.io/github/carlosp420/seqrecord-expanded?branch=master

.. image:: https://www.quantifiedcode.com/api/v1/project/b0bf8d6e31704c11abeb0b9043c11891/badge.svg
   :target: https://www.quantifiedcode.com/app/project/b0bf8d6e31704c11abeb0b9043c11891
   :alt: Code issues

seqrecord-expanded
==================
Another SeqRecord class with methods: degenerate seqs, codon positions based on
reading frames, etc.

Usage
-----
By default it assumes a DNA sequence with ambiguous characters.

.. code:: python

    >>> from seqrecord_expanded import SeqRecordExpanded
    >>> seq_record = SeqRecordExpanded('TCTGAATGGAAGACAAAGCGTCCA',
    ...                                voucher_code='CP100-09',
    ...                                taxonomy={'genus': 'Melitaea',
    ...                                          'species': 'phoebe',
    ...                                         },
    ...                                gene_code='EF1a',
    ...                                reading_frame=1,
    ...                                table=1,  # translation table
    ...                                )
    >>> # Degenerate sequence standard genetic code
    >>> seq_record.degenerate()
    'TCNGARTGGAARACNAARMGNCCN'
    >>>
    >>> # Degenerate sequence S method
    >>> seq_record.degenerate(method='S')
    'AGYGARTGGAARACNAARMGNCCN'
    >>>
    >>> # Degenerate sequence Z method
    >>> seq_record.degenerate(method='Z')
    'TCNGARTGGAARACNAARMGNCCN'
    >>>
    >>> # Degenerate sequence SZ method
    >>> seq_record.degenerate(method='SZ')
    'NNNGARTGGAARACNAARMGNCCN'
    >>>
    >>> # get first codon positions
    >>> seq_record.first_codon_position()
    'TGTAAACC'
    >>>
    >>> # get second codon positions
    >>> seq_record.second_codon_position()
    'CAGACAGC'
    >>>
    >>> # get third codon positions
    >>> seq_record.third_codon_position()
    'TAGGAGTA'
    >>>
    >>> # get first and second positions
    >>> seq_record.first_and_second_positions()
    'TCGATGAAACAACGCC'
    >>>
    >>> # translate
    >>> seq_record.translate()
    OKPDOR

Installation
------------

.. code-block:: shell

    pip install seqrecord-expanded

Requirements
^^^^^^^^^^^^
Python3.

.. code-block:: shell

    pip install -r requirements.txt


Compatibility
-------------

Licence
-------
BSD.

Authors
-------

`seqrecord-expanded` was written by `Carlos Peña <mycalesis@gmail.com>`_.
