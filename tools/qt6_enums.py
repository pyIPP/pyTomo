"""Rewrite short Qt5 enum access (Qt.AlignCenter, QSizePolicy.Expanding) to the fully
scoped form (Qt.AlignmentFlag.AlignCenter) accepted by PyQt6 and PyQt5 >= 5.11.

usage: python tools/qt6_enums.py FILE [FILE ...]   (needs PyQt6; edits in place)
"""
import ast
import enum
import sys

from PyQt6 import QtCore, QtGui, QtWidgets

MODULES = (QtCore, QtGui, QtWidgets)


def scoped(cls_name, member):
    for m in MODULES:
        cls = getattr(m, cls_name, None)
        if cls is None or hasattr(cls, member):
            continue
        for name in dir(cls):
            e = getattr(cls, name)
            if isinstance(e, type) and issubclass(e, enum.Enum) and member in e.__members__:
                return '%s.%s.%s' % (cls_name, name, member)
    return None


def convert(path):
    src = open(path).read()
    edits, missing = [], set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) \
                and node.value.id.startswith('Q'):
            new = scoped(node.value.id, node.attr)
            if new:
                edits.append((node.lineno, node.col_offset, node.end_col_offset, new))
    lines = src.splitlines(keepends=True)
    for ln, c0, c1, new in sorted(edits, reverse=True):
        lines[ln - 1] = lines[ln - 1][:c0] + new + lines[ln - 1][c1:]
    open(path, 'w').write(''.join(lines))
    return len(edits)


if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(p, convert(p), 'enums scoped')
