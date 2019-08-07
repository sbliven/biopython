"""Common structures for parsing MolQL syntaxes

"""


class MolQLFormatError(ValueError):
    """Indicates an error parsing a MolQL query"""
    pass


class Expression(object):
    """Represents an unevaluated MolQL expression

    ```
    Expression :: Literal | Apply
    Literal :: string | number | boolean
    Apply :: Head Arguments?
    Head :: Symbol
    Arguments :: Expression+
    ```

    Since literal values are parsed directly to their python equivalents, an
    Expression object always represents the 'Apply' form and contains a *head*
    (the function name) and zero or more *args* (the argument)

    """
    def __init__(self, head, args=[]):
        self.head = head
        self.args = args

    def __repr__(self):
        return 'Expression("%s"%s)' % (self.head, ", ..." if len(self.args) > 0 else "")
