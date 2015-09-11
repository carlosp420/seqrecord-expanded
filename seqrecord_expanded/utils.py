import itertools


def chain_and_flatten(seq1, seq2):
    """Takes two strings (first and second codon positions) and chains and
    intercalates them.

    :returns: string of combining the two seq strings.
    """
    my_chain = itertools.zip_longest(seq1, seq2)
    out = [i for i in itertools.chain.from_iterable(my_chain) if i]
    return ''.join(out)
