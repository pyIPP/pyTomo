"""Replace star imports of numpy / matplotlib.pyplot / pylab by explicit module aliases.

    from numpy import *            ->  import numpy as np        (+ np.name)
    from matplotlib.pyplot import * ->  import matplotlib.pyplot as plt  (+ plt.name)
    from pylab import *            ->  both of the above (+ np.linalg., np.random., ...)

Semantics are preserved for the reference environment (NumPy 1.26, Matplotlib 3.8):
tools/star_names.json (made by record_star_names.py) maps each exported name to its
origin, and builtins (max, min, abs, round in pylab) stay builtins. NumPy 2 exports
max/min/abs/round/pow/bool as well, which silently changed the meaning of
`from numpy import *`.

Names bound in the file (assignment, parameter, import, `except ... as e`) are not
rewritten and are reported, because the star name might be used in another scope.
Check the result with pyflakes and the test suite.

usage: python tools/dewildcard.py FILE [FILE ...]   (edits files in place)
"""
import ast
import json
import os
import sys
import symtable

HERE = os.path.dirname(os.path.abspath(__file__))
STAR = json.load(open(os.path.join(HERE, 'star_names.json')))
IMPORTS = {'np': 'import numpy as np', 'plt': 'import matplotlib.pyplot as plt',
           'mlab': 'from matplotlib import mlab'}


def _bound_names(table, out):
    for sym in table.get_symbols():
        if sym.is_assigned() or sym.is_imported() or sym.is_parameter():
            out.add(sym.get_name())
    for child in table.get_children():
        _bound_names(child, out)
    return out


def convert(path):
    src = open(path).read()
    tree = ast.parse(src)
    stars = [n for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
             and n.module in STAR and any(a.name == '*' for a in n.names)]
    if not stars:
        return 0, [], []
    names = {}
    for node in sorted(stars, key=lambda n: n.lineno):  # later star imports win
        names.update(STAR[node.module])
    bound = _bound_names(symtable.symtable(src, path, 'exec'), set())

    edits, skipped, unknown, used = [], set(), set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id in names:
            qual = names[node.id]
            if node.id in bound:
                skipped.add(node.id)
            elif qual is None:
                unknown.add(node.id)
            elif qual != node.id:
                edits.append((node.lineno, node.col_offset, len(node.id), qual))
                used.add(qual.split('.')[0])
    lines = src.splitlines(keepends=True)
    for lineno, col, n, qual in sorted(edits, reverse=True):
        line = lines[lineno - 1]
        lines[lineno - 1] = line[:col] + qual + line[col + n:]

    # replace the star import lines by the imports actually needed
    top = {ast.unparse(n) for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))}
    first = True
    for node in sorted(stars, key=lambda n: n.lineno):
        line = lines[node.lineno - 1]
        indent = line[:len(line) - len(line.lstrip())]
        new = ''
        if first:
            new = ''.join(indent + IMPORTS[a] + '\n' for a in sorted(used)
                          if not (indent == '' and IMPORTS[a] in top))
            first = False
        lines[node.lineno - 1] = new
    open(path, 'w').write(''.join(lines))
    return len(edits), sorted(skipped), sorted(unknown)


if __name__ == '__main__':
    for p in sys.argv[1:]:
        n, skipped, unknown = convert(p)
        msg = '%s: %d names qualified' % (p, n)
        if skipped:
            msg += '\n   review (also bound locally): ' + ', '.join(skipped)
        if unknown:
            msg += '\n   unknown origin, left as is: ' + ', '.join(unknown)
        print(msg)
