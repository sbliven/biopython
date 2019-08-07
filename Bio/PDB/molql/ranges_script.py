"""Parse 'ranges' MolQL format

"""


try:
    import pyparsing as pp
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError("Install pyparsing to use MolQL "
                                       "(e.g. pip install pyparsing)")


# Singleton
_ranges_bnf = None


def rangesBNF():
    """Construct the grammar for the "range" format

    Returns a grammar, which upon successfully parsing should produce a list of ResRange.

    """
    global _ranges_bnf
    if _ranges_bnf is None:
        def parseResi(s, l, t):
            if len(t) == 1:
                return (t[0], " ")
            else:
                return tuple(t)

        def parseRange(s, l, t):
            return ResRange(*t)

        inscode = pp.Word(pp.alphanums)
        signedInt = pp.Regex("[+-]?[0-9]+").setParseAction(lambda s, l, t: [int(t[0])])
        resi = (signedInt + pp.Optional(inscode)).leaveWhitespace().setParseAction(parseResi)
        chain_id = pp.Word(pp.alphanums)
        rangeExpr = (chain_id.setResultsName('chain_id') +
                     pp.Optional(pp.Suppress('.') + resi.setResultsName('start') +
                                 pp.Optional(pp.Suppress('-') + resi.setResultsName('end')))
                     ).leaveWhitespace().setParseAction(parseRange)
        ranges = pp.delimitedList(rangeExpr)
        _ranges_bnf = ranges
    return _ranges_bnf


class ResRange(object):
    """Represents a simple residue range

    Ranges are inclusive between start and end. If one endpoint is missing
    than the beginning or end of the chain is implied

    Arguments:
    - chain_id (str): Author chain ID
    - start (int or (int,str)): starting residue number, optionally with insertion code
    - end (int or (int,str)): ending residue number, optionally with insertion code
    """
    def __init__(self, chain_id, start=None, end=None):
        self.chain_id = chain_id
        if type(start) is tuple:
            self.start = start
        else:
            self.start = (start, " ")
        if type(end) is tuple:
            self.end = end
        else:
            self.end = (end, " ")

    def __str__(self):
        return "%s.%s%s-%s%s" % (self.chain_id,
                                 "^" if self.start is None or self.start[0] is None else self.start[0],
                                 "" if self.start is None or self.start[1] is None else self.start[1].strip(),
                                 "$" if self.end is None or self.end[0] is None else self.end[0],
                                 "" if self.end is None or self.end[1] is None else self.end[1].strip())

    def __repr__(self):
        return str(self)
