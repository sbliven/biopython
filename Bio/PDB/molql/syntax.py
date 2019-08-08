"""Shared structures for reading and writing MolQL in different syntaxes

To add a syntax, subclass Syntax with load/dump methods. Then call 'register_syntax'
once in your module.
"""
import abc
from StringIO import StringIO

# Registered syntaxes. {name (str): Syntax}
_syntaxes = {}


def register_syntax(syntax):
    """Register a new MolQL syntax.

    `syntax.name()` should be unique. Re-registering a syntax with the same name
    will cause that syntax to be replaced. Names are case-insensitive

    Arguments:
    - syntax (Syntax): singleton for the new syntax
    """
    _syntaxes[syntax.name().lower()] = syntax


def get_syntax(name):
    """Get a syntax by name

    Arguments:
    - name (str): syntax name (case insensitive)

    Returns: (Syntax)
    """
    return _syntaxes[name.lower()]


def get_syntaxes():
    """Get a list of all registered syntaxes

    Returns: sequence of str
    """
    return _syntaxes.keys()


class Syntax(object):
    """Abstract class representing a MolQL syntax

    Individual function calls should be independent, so the class should generally
    not maintain internal state.

    A singleton for new syntaxes should be registered using `register_syntax`.

    """
    __metaclass__ = abc.ABCMeta

    @staticmethod
    @abc.abstractmethod
    def name(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def load(self, fp):
        """Parse query from a file-like object

        Arguments:
        - fp (file-like): query according to this syntax

        Returns: Expression
        """
        raise NotImplementedError()

    def loads(self, s):
        """Parse query from a string

        Arguments:
        - s (str): query according to this syntax

        Returns: Expression
        """
        fp = StringIO(s)
        return self.load(fp)

    def dump(self, query, fp):
        """Write a MolQL query to a file

        This method is optional. Syntaxes that do not support writing should
        raise NotImplementedError.

        If the query contains operations not supported by this syntax, raise
        MolQLFormatError.

        Arguments:
        - query (Expression): query to write
        - fp (file-like): output file

        """
        raise NotImplementedError("Writing is not supported by this syntax")

    def dumps(self, query):
        """Convert a MolQL query to a string following this syntax

        This method is optional. Syntaxes that do not support writing should
        raise NotImplementedError.

        If the query contains operations not supported by this syntax, raise
        MolQLFormatError.

        Arguments:
        - query (Expression): query to write

        Returns: (str)

        """
        s = StringIO(query)
        self.dump(query, s)
        return s.getvalue()


class MolQLFormatError(ValueError):
    """Indicates an error parsing a MolQL query"""
    pass


class Expression(object):
    """Represents an unevaluated MolQL expression

    ```
    Expression :: Literal | Apply
    Literal :: string | number | boolean
    Apply :: Head Argument*
    Head :: Symbol
    Argument :: (Symbol ':')? Expression
    ```

    Since literal values are parsed directly to their python equivalents, an
    Expression object always represents the 'Apply' form and contains a *head*
    (the function name) and zero or more *args* (the argument)

    """
    def __init__(self, head, args={}):
        """
        Arguments:
        - head (str): name of operation
        - args ({str: Expression}): list of arguments. Positional arguments
          are represented by integer keys; named arguments have arbitrary keys.
        """
        self.head = head
        self.args = args

    def __repr__(self):
        return 'Expression("%s"%s)' % (self.head, ", ..." if len(self.args) > 0 else "")
