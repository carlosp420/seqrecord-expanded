import setuptools

setuptools.setup(
    name="seqrecord-expanded",
    version="0.0.0",
    url="https://github.com/carlosp420/seqrecord-expanded",

    author="Carlos Peña",
    author_email="mycalesis@gmail.com",

    description="BioPython's SeqRecord class, but expanded with additional methods: degenerate seqs, codon positions based on reading frames, etc.",
    long_description=open('README.rst').read(),

    packages=['seqrecord_expanded'],

    install_requires=[],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],
)
