seqrecord-expanded
==================

.. list-table::
    :stub-columns: 1

    * - tests
      - | |travis| |requires| |coveralls|
        | |quantified-code|
    * - package
      - |version| |wheel| |supported-versions| |supported-implementations|

.. |travis| image:: https://travis-ci.org/carlosp420/seqrecord-expanded.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/carlosp420/seqrecord-expanded

.. |requires| image:: https://requires.io/github/carlosp420/seqrecord-expanded/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/carlosp420/seqrecord-expanded/requirements/?branch=master

.. |coveralls| image:: https://coveralls.io/repos/carlosp420/seqrecord-expanded/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/carlosp420/seqrecord-expanded

.. |version| image:: https://img.shields.io/pypi/v/seqrecord-expanded.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/seqrecord-expanded

.. |wheel| image:: https://img.shields.io/pypi/wheel/seqrecord-expanded.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/seqrecord-expanded

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/seqrecord-expanded.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/seqrecord-expanded

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/seqrecord-expanded.svg?style=flat
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/seqrecord-expanded

.. |quantified-code| image:: https://www.quantifiedcode.com/api/v1/project/b0bf8d6e31704c11abeb0b9043c11891/badge.svg
   :alt: Code issues
   :target: https://www.quantifiedcode.com/app/project/b0bf8d6e31704c11abeb0b9043c11891


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
    >>> # translate to aminoacid sequence
    >>> seq_record.translate()
    'SEWKTKRP'
    >>> # translate to aminoacid sequence
    >>> seq_record.translate(table=1)
    'SEWKTKRP'

Installation
------------

.. code-block:: shell

    pip install seqrecord-expanded

Requirements
^^^^^^^^^^^^

.. code-block:: shell

    pip install -r requirements.txt


Compatibility
-------------
Supported Python versions: 2.6, 2.7, 3.3, 3.4, 3.5, pypy.

Licence
-------
BSD.

Authors
-------

`seqrecord-expanded` was written by `Carlos Peña <mycalesis@gmail.com>`_.
