import itertools

import six
if six.PY2:
    from six.moves import zip_longest
else:
    from itertools import zip_longest


def chain_and_flatten(seq1, seq2):
    """Takes two strings (first and second codon positions) and chains and
    intercalates them.

    Returns:
        (str): String of combining the two seq strings.
    """
    my_chain = zip_longest(seq1, seq2)
    out = [i for i in itertools.chain.from_iterable(my_chain) if i]
    return ''.join(out)
