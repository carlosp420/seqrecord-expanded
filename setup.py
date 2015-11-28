#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import io
import setuptools
from os.path import join
from os.path import dirname


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


setuptools.setup(
    name="seqrecord_expanded",
    version="0.2.3",
    license="BSD",
    url="https://github.com/carlosp420/seqrecord-expanded",

    author="Carlos Peña",
    author_email="mycalesis@gmail.com",

    description="Another SeqRecord class with methods: degenerate seqs, codon positions based on reading frames, etc.",
    long_description='%s' % read('README.rst'),

    packages=['seqrecord_expanded'],

    install_requires=[
        'biopython==1.66',
        'degenerate-dna==0.0.9',
        'six==1.10.0',
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
