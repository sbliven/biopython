"""Parse MolQL's lisp-like format

"""
from . import syntax


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
        # TODO Stub
        raise NotImplementedError()

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
