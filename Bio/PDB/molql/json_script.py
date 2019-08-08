"""Parse MolQL queries from JSON syntax

"""
import json
from ... import __version__
from . import syntax


class JSONSyntax(syntax.Syntax):
    @staticmethod
    def name():
        return "json"

    def loads(self, query):
        return self.load_dict(json.loads(query))

    def load(self, fp):
        return self.load_dict(json.load(fp))

    def load_dict(self, jsonDict):
        """Parse a MolQL query from JSON dictionary.

        Arguments:
        - jsonDict (dict): raw dictionary representation of the JSON

        Returns: Expression
        """
        # Currently ignore source and version information
        if "expression" not in jsonDict:
            raise syntax.MolQLFormatError("Illegal MolQL query (json). Missing 'expression'")

        def parse_expression(node):
            if isinstance(node, dict):
                if "head" not in node:
                    raise syntax.MolQLFormatError("Illegal MolQL query (json). Missing 'head'")
                if "args" not in node:
                    return syntax.Expression(node["head"])
                else:
                    if isinstance(node["args"], dict):
                        return syntax.Expression(node["head"],
                                                 {name: parse_expression(n)
                                                  for name, n in node["args"].items()})
                    elif isinstance(node["args"], list):
                        return syntax.Expression(node["head"],
                                                 {str(i): parse_expression(n)
                                                  for i, n in enumerate(node["args"])})
                    else:
                        raise ValueError("Error Parsing MolQL (JSON): Unexpected args value")

            else:
                return node  # leaf

        return parse_expression(jsonDict["expression"])

    def dump_dict(self, expr):
        """Convert Query to JSON

        Arguments:
        - molql (Expression): query
        Returns: (str)
        """
        def to_json_dict(tree):
            d = {"head": tree.head}
            if len(tree.args) > 0:
                d["args"] = {label: to_json_dict(child) if isinstance(child, syntax.Expression) else child
                             for label, child in tree.args.items()}
            return d

        d = {"source": "BioPython %s" % __version__,
             "version": "0.1.0",
             "expression": to_json_dict(expr)
             }
        return d

    def dumps(self, expr, **kw):
        return json.dumps(self.dump_dict(expr), **kw)

    def dump(self, expr, fp, **kw):
        json.dump(self.dump_dict(expr), fp, **kw)


syntax.register_syntax(JSONSyntax())
