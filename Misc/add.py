#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

# Ensure arith.py is importable from ..\Arith relative to this script
arith_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "Arith")
)
sys.path.insert(0, arith_dir)

import arith


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  add.py input.fit output.fit operand [offset]\n"
        "\n"
        "Equivalent to:\n"
        "  arith.py add input.fit output.fit operand [offset]\n"
    )
    sys.exit(1)


def main():
    if len(sys.argv) not in (4, 5):
        usage()

    # Rewrite argv for arith.py:
    # arith.py add input output operand [offset]
    new_argv = [sys.argv[0], "add"] + sys.argv[1:]
    sys.argv = new_argv
    arith.main()


if __name__ == "__main__":
    main()
