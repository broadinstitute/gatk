#!/usr/bin/env python3
"""Re-pad GitHub-flavored Markdown tables so every row is rectangular.

IntelliJ's Markdown plugin reformats tables to rectangular form on save, so a table
committed with ragged rows produces spurious diffs the next time anyone opens the file
in the IDE. A single cell one character wider than its column is enough to do it, and
that is nearly invisible in review. Running this after editing a table keeps the
on-disk form stable no matter who or what edited it last.

Formatting matches what the IDE produces: each cell is padded to the widest content in
its column, with one space of padding on each side, and separator rows are filled with
dashes across the full cell width. Alignment markers (`:---`, `---:`, `:---:`) are
preserved. Content inside fenced code blocks is never touched, so embedded JSON,
SQL and shell examples are safe.

Usage:
    format_markdown_tables.py FILE_OR_DIR [FILE_OR_DIR ...]   # rewrite in place
    format_markdown_tables.py --check FILE_OR_DIR [...]       # report only, exit 1 if unformatted

Directories are searched recursively for *.md. Exit status is 0 when nothing needed
changing (or everything was rewritten), 1 in --check mode when a file would change.
"""

import argparse
import os
import re
import sys

FENCE = re.compile(r'^\s*(```|~~~)')
SEPARATOR_CELL = re.compile(r'^:?-+:?$')
UNESCAPED_PIPE = re.compile(r'(?<!\\)\|')


def _split_cells(line):
    """Cells of a pipe table row.

    Drops the empty piece before the leading pipe, and the one after the trailing pipe
    only when the row actually has a trailing pipe — a row written without one still
    carries content in its final cell, which must not be discarded. Escaped pipes
    (``\\|``) do not separate cells.
    """
    stripped = line.strip()
    pieces = UNESCAPED_PIPE.split(stripped)[1:]
    if stripped.endswith('|') and not stripped.endswith('\\|') and pieces:
        pieces = pieces[:-1]
    return [piece.strip() for piece in pieces]


def _is_separator_row(line):
    cells = _split_cells(line)
    return bool(cells) and all(SEPARATOR_CELL.match(cell) for cell in cells)


def _separator_cell(width, marker):
    """Render a separator cell occupying exactly ``width + 2`` characters.

    Alignment colons are preserved and counted against the width, so the separator row
    stays flush with the content rows; padding a separator to a minimum dash count
    instead would reintroduce the raggedness this script exists to remove.
    """
    left = marker.startswith(':')
    right = marker.endswith(':') and marker != ':'
    dashes = width + 2 - int(left) - int(right)
    return (':' if left else '') + '-' * max(dashes, 1) + (':' if right else '')


def format_table(rows):
    """Given the lines of one pipe table, return them padded rectangularly."""
    parsed = [_split_cells(row) for row in rows]
    ncols = max(len(cells) for cells in parsed)
    parsed = [cells + [''] * (ncols - len(cells)) for cells in parsed]

    # Minimum width of 1 keeps an empty column's separator to a legal single dash.
    widths = [1] * ncols
    for index, cells in enumerate(parsed):
        if index == 1:
            continue  # the separator row imposes no content width
        for column, cell in enumerate(cells):
            widths[column] = max(widths[column], len(cell))

    out = []
    for index, cells in enumerate(parsed):
        if index == 1:
            out.append('|' + '|'.join(
                _separator_cell(width, marker) for width, marker in zip(widths, cells)
            ) + '|')
        else:
            out.append('|' + '|'.join(
                ' ' + cell.ljust(width) + ' ' for cell, width in zip(cells, widths)
            ) + '|')
    return out


def format_text(text):
    """Return (formatted_text, tables_changed)."""
    lines = text.split('\n')
    out, index, in_fence, changed = [], 0, False, 0

    while index < len(lines):
        line = lines[index]

        if FENCE.match(line):
            in_fence = not in_fence
            out.append(line)
            index += 1
            continue

        if not in_fence and line.startswith('|'):
            end = index
            while end < len(lines) and lines[end].startswith('|'):
                end += 1
            table = lines[index:end]
            # A pipe table needs a header and a separator row; anything else is left alone.
            if len(table) >= 2 and _is_separator_row(table[1]):
                formatted = format_table(table)
                if formatted != table:
                    changed += 1
                out.extend(formatted)
                index = end
                continue

        out.append(line)
        index += 1

    return '\n'.join(out), changed


def markdown_files(paths):
    for path in paths:
        if os.path.isdir(path):
            for root, _, names in os.walk(path):
                for name in sorted(names):
                    if name.endswith('.md'):
                        yield os.path.join(root, name)
        elif path.endswith('.md'):
            yield path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('paths', nargs='+', help='Markdown files, or directories to search recursively')
    parser.add_argument('--check', action='store_true',
                        help='report files that would change without rewriting them')
    parser.add_argument('--quiet', action='store_true', help='only report problems')
    args = parser.parse_args(argv)

    would_change = []
    for path in markdown_files(args.paths):
        with open(path) as handle:
            original = handle.read()
        formatted, changed = format_text(original)
        if formatted == original:
            continue
        would_change.append(path)
        if args.check:
            print(f'{path}: {changed} table(s) not rectangular')
        else:
            with open(path, 'w') as handle:
                handle.write(formatted)
            if not args.quiet:
                print(f'{path}: {changed} table(s) reformatted')

    if args.check and would_change:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
