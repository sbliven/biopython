"""Parse MolQL's lisp-like format

"""
try:
    import pyparsing as pp
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError("Install pyparsing to use molql-script "
                                       "(e.g. pip install pyparsing)")

from . import syntax


# Singleton
_bnf = None


def molscriptBNF():
    """Construct the grammar for the "molql" format

    Returns a grammar, which upon successfully parsing should produce an Expression.

    """
    global _bnf
    if _bnf is None:

        def parseApply(s, l, t):
            print("Creating Expression(%s)" % (t))
            if len(t) == 1:
                return syntax.Expression(t[0])
            else:
                return syntax.Expression(t[0], t[1])

        def parseArgument(s, l, t):
            "Parse Argument as (label, Expression)"
            print("Creating Argument(%s)" % (t))
            if len(t) == 1:
                return (None, t[0])
            else:
                return (t[0], t[1])

        def parseArguments(s, l, t):
            args = dict()
            for i, pair in enumerate(t):
                label, expr = pair
                args[i if label is None else label] = expr
            print("Creating Arguments(%s)" % args)
            return args

        # _e expression
        # _s symbol
        token = pp.Regex(r"[^\s(){}'\":]+").setParseAction(lambda s, l, t: t[0])
        string_s = pp.quotedString.setParseAction(pp.removeQuotes)
        # Allow bare strings as fallback. Validate literals after parsing
        bare_string_s = token.copy()
        fnumber_s = pp.Regex(r"[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?").setParseAction(lambda s, l, t: float(t[0]))
        bool_s = (pp.CaselessKeyword("true").setParseAction(lambda s, l, t: True) |
                  pp.CaselessKeyword("false").setParseAction(lambda s, l, t: False))

        literal_e = pp.Group(string_s | fnumber_s | bool_s | bare_string_s)

        # Use non-back-tracing '-' operator for better error handling
        arg_name_s = (pp.Suppress(":") - token)
        expr_e = pp.Forward()
        kwargument_e = (arg_name_s - expr_e).setParseAction(parseArgument)
        # How do we get pyparsing to treat this as a new expression?
        posargument_e = expr_e.copy().setParseAction(parseArgument)
        arguments_e = (pp.ZeroOrMore(posargument_e) -
                       pp.ZeroOrMore(kwargument_e)).setParseAction(parseArguments)

        def headAction(s, l, t):
            print("Head: %s" % (t))
            return t

        head_s = token.copy().setParseAction(headAction)
        apply_e = (pp.Suppress("(") -
                   head_s -
                   arguments_e -
                   pp.Suppress(")")
                   ).setParseAction(parseApply)
        expr_e <<= (apply_e | literal_e).setParseAction(lambda s, l, t: t[0])

        _bnf = expr_e
    return _bnf


class MolQLSyntax(syntax.Syntax):
    @staticmethod
    def name():
        return "molql"

    def load(self, fp):
        """Parse query from a file-like object

        Arguments:
        - fp (file-like): query according to this syntax

        Returns: Expression
        """
        return molscriptBNF().parseFile(fp, True)[0]

    def loads(self, s):
        """Parse query from a string

        Arguments:
        - s (str): query according to this syntax

        Returns: Expression
        """
        return molscriptBNF().parseString(s)[0]

    def dump(self, query, fp, indent=None):
        """Write a MolQL query to a file using lisp-like syntax

        Supports all MolQL operations

        Arguments:
        - query (Expression): query to write
        - fp (file-like): output file
        - indent (int?): If not None, pretty print the output with the initial
          indentation specified. Arguments will be indented by 2 additional spaces.

        """

        pretty = indent is not None
        if pretty:
            fp.write(" " * indent)
        fp.write("(")
        fp.write(query.head)
        if len(query.args) > 0:
            if pretty:
                fp.write("\n")
            else:
                fp.write(" ")
            for label, child in query.args.items():
                if isinstance(child, syntax.Expression):
                    self.dumps(child, fp, indent=indent + 2)
                else:
                    if pretty:
                        fp.write(" " * (indent + 2))
                        fp.write(repr(child))
                        fp.write("\n")
            if pretty:
                fp.write(" " * indent)
        fp.write(")")
        if pretty:
            fp.write("\n")


syntax.register_syntax(MolQLSyntax())
