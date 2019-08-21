"""Registry of implementations for MolQL symbols.

This maps a symbol name (e.g. "structure.generator.atom-groups") to an implementing function.
"""
# All known language primitives
_symbols = {}


def register_symbols(symbols):
    """Register a set of MolQL symbols.

    Symbols implement particular MolQL keywords

    Args:
    - symbols (Dict[str, Callable]): Map from the symbol name to a Callable
      implementing that symbol.
    """
    _symbols.update(symbols)


def get_symbol(symbol_name):
    """Look up a symbol by name"""
    return _symbols[symbol_name]


def symbol(symbol_name):
    """Decorator to declare a MolQL symbol

    The following are equivalent:

        @symbol("=") def eq(a, b): a == b

    and

        def eq(a, b): a == b
        register_symbols({"=": eq})

    """
    def symbol_decorator(fn):
        _symbols[symbol_name] = fn
        return fn
    return symbol_decorator