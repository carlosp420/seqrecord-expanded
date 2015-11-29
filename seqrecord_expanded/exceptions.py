class MissingParameterError(Exception):
    pass


class TranslationErrorMixedGappedSeq(Exception):
    def __init__(self, voucher_code, gene_code, e):
        self.voucher_code = voucher_code
        self.gene_code = gene_code
        self.msg = e.__str__().replace("N", "?")

    def __str__(self):
        return "Translation Error in sequence {0!r} {1!r}: {2}.".format(self.voucher_code,
                                                                        self.gene_code,
                                                                        self.msg)
