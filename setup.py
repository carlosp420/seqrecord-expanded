# -*- encoding: utf-8 -*-
import setuptools


setuptools.setup(
    name="seqrecord-expanded",
    version="0.1.1",
    license="BSD",
    url="https://github.com/carlosp420/seqrecord-expanded",

    author="Carlos Peña",
    author_email="mycalesis@gmail.com",

    description="Another SeqRecord class with methods: degenerate seqs, codon positions based on reading frames, etc.",
    long_description=open('README.rst').read(),

    packages=['seqrecord_expanded'],

    install_requires=[
        'biopython==1.65',
        'degenerate-dna==0.0.9',
        'six==1.9.0',
    ],

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='tests',
)
