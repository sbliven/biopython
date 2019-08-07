"""Parse MolQL's lisp-like format

"""

from .parsing import Expression


def to_molql(self, indent=None):
    """Output in MolQL lisp-like format"""
    pretty = indent is not None
    tokens = []
    if pretty:
        tokens.append(" " * indent)
    tokens.extend(("(", self.head))
    if len(self.children) > 0:
        if pretty:
            tokens.append("\n")
        else:
            tokens.append(" ")
        for child in self.children:
            if isinstance(child, Expression):
                tokens.extend(child.to_molql(pretty, indent + 2))
            else:
                if pretty:
                    tokens.append(" " * (indent + 2))
                    tokens.append(repr(child))
                    tokens.append("\n")
        if pretty:
            tokens.append(" " * indent)
    tokens.append(")")
    if pretty:
        tokens.append("\n")

    return "".join(tokens)
