"""Parse MolQL queries from JSON syntax

"""

import json
from .parsing import Expression, MolQLFormatError
from ... import __version__


def parse_json(query):
    """Parse a MolQL query from JSON format.

    Arguments:
    - json (string or file-like): json input

    Returns: MolQL
    """
    if isinstance(query, str):
        jsonDict = json.loads(query)
    else:
        jsonDict = json.load(query)

    # Currently ignore source and version information
    if "expression" not in jsonDict:
        raise MolQLFormatError("Illegal MolQL query (JSON). Missing 'expression'")

    def parse_expression(node):
        if isinstance(node, dict):
            if "head" not in node:
                raise MolQLFormatError("Illegal MolQL query (JSON). Missing 'head'")
            if "args" not in node:
                return Expression(node["head"])
            else:
                return Expression(node["head"], [parse_expression(n) for n in node["args"]])
        else:
            return node  # leaf

    return parse_expression(jsonDict["expression"])


def to_json(molql, indent=None):
    """Convert Query to JSON

    Arguments:
    - molql (MolQL): query
    Returns: (str)
    """
    def to_json_dict(tree):
        d = {"head": tree.head}
        if len(tree.children) > 0:
            d["args"] = [child._to_json_dict() if isinstance(child, Expression) else child
                         for child in tree.children]
        return d

    d = {"source": "BioPython %s" % __version__,
         "version": "0.1.0",
         "expression": to_json_dict(molql)
         }
    return json.dumps(d, indent=indent)
