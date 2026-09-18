#!/usr/bin/env python3
"""Unit tests for the check-wdl-heredocs helper.

The failure this guards against is expensive to find any other way: `womtool validate`
passes, `bash -n` passes for half of it, and the break surfaces only when the task runs in
the cloud, reported at the generated script's last line rather than anywhere near its
cause. So the fixtures below are written as the two failure modes actually appear in a
diff, and both directions matter -- a checker that reports correct code gets switched off,
and this repo contains WDLs written under both of the valid conventions.

The helper is deliberately extensionless -- see its docstring for why -- so it is loaded by
path rather than imported by name.
"""

import contextlib
import importlib.machinery
import importlib.util
import io
import pathlib
import tempfile
import unittest

_HELPER = pathlib.Path(__file__).resolve().parents[1] / 'check-wdl-heredocs'
_spec = importlib.util.spec_from_loader(
    'check_wdl_heredocs',
    importlib.machinery.SourceFileLoader('check_wdl_heredocs', str(_HELPER)),
)
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)


def wdl(block):
    """A minimal WDL whose single command block holds `block` verbatim.

    The lines are not re-indented, because their indentation is the entire subject here.
    """
    return '\n'.join(['task T {', '    command <<<', *block, '    >>>', '}', ''])


def broken(text):
    """Only the BROKEN findings.

    A fixture that has to defeat the dedent needs a line at column zero, which is itself
    reported FRAGILE. That is correct -- such a block really is never dedented -- so tests
    about the breakage filter for it rather than asserting on the whole list and quietly
    encoding the incidental finding as expected.
    """
    return [f for f in checker.findings(text) if f.severity == 'BROKEN']


def fragile(text):
    return [f for f in checker.findings(text) if f.severity == 'FRAGILE']


class TestValidConventions(unittest.TestCase):
    """Neither convention may be reported broken. Both are in use in this repo."""

    def test_terminator_indented_with_the_block_is_clean(self):
        """The dedent strips the common prefix from body and terminator alike."""
        self.assertEqual([], checker.findings(wdl([
            '        echo hello',
            "        python3 <<'CODE'",
            '        import sys',
            '        print(sys.argv)',
            '        CODE',
        ])))

    def test_everything_at_column_zero_is_fragile_but_not_broken(self):
        """The other convention: no dedent happens, so the source is already correct."""
        found = checker.findings(wdl([
            'echo hello',
            "python3 <<'CODE'",
            'import sys',
            'print(sys.argv)',
            'CODE',
        ]))
        self.assertEqual(['FRAGILE'], [f.severity for f in found])
        self.assertIn('never dedented', found[0].message)

    def test_a_body_that_needs_its_indentation_is_left_alone(self):
        """Nested Python has to keep its own structure, and the dedent is uniform."""
        self.assertEqual([], checker.findings(wdl([
            "        python3 <<'CODE'",
            '        if True:',
            '            print("nested")',
            '        CODE',
        ])))


class TestTerminatorInvariant(unittest.TestCase):
    """A terminator not at column zero after the dedent: bash reads to end of file."""

    def test_one_column_zero_line_defeats_the_dedent(self):
        found = checker.findings(wdl([
            '        echo hello',
            'this_line_is_at_column_zero=1',
            "        python3 <<'CODE'",
            '        import sys',
            '        CODE',
        ]))
        self.assertEqual(['BROKEN'], [f.severity for f in found])
        self.assertIn('CODE', found[0].message)
        self.assertIn('no terminator at column zero', found[0].message)

    def test_every_affected_heredoc_is_named_not_just_the_first(self):
        """bash reports one. Fixing only that one leaves the next still broken."""
        found = checker.findings(wdl([
            'at_column_zero=1',
            "        python3 <<'FIRST'",
            '        import sys',
            '        FIRST',
            '        cat > x.sql <<SECOND',
            '        select 1;',
            '        SECOND',
        ]))
        self.assertEqual(1, len(found))
        self.assertIn('FIRST', found[0].message)
        self.assertIn('SECOND', found[0].message)

    def test_dash_heredocs_are_exempt(self):
        """`<<-` strips leading tabs, so its terminator may be indented."""
        self.assertEqual([], broken(wdl([
            'at_column_zero=1',
            '        cat <<-CODE',
            '        indented body',
            '        CODE',
        ])))


class TestBodyInvariant(unittest.TestCase):
    """The quieter half: the terminator is where bash wants it, the body is not.

    This is what the obvious repair for a terminator error produces -- outdent the
    terminator, leave the body alone -- and it trades `unexpected end of file` for
    `IndentationError` at `"<stdin>", line 1`, naming neither the WDL nor the task.
    """

    HALF_FIXED = [
        'at_column_zero=1',
        "        python3 <<'CODE'",
        '        import sys',
        '        print(sys.argv)',
        'CODE',
    ]

    def test_an_indented_python_body_is_broken(self):
        found = broken(wdl(self.HALF_FIXED))
        self.assertEqual(1, len(found))
        self.assertIn('still indented by 8', found[0].message)

    def test_it_points_at_the_line_that_opens_the_heredoc(self):
        """Not at the block, and not at the IndentationError's own line 1.

        `task T {` is line 1 and `command <<<` line 2, so HALF_FIXED's opener is line 4.
        """
        self.assertEqual(4, broken(wdl(self.HALF_FIXED))[0].line)

    def test_the_column_zero_line_that_caused_it_is_reported_too(self):
        """Both are true and both are worth printing: one is the breakage, one the cause."""
        self.assertEqual(1, len(fragile(wdl(self.HALF_FIXED))))

    def test_an_indented_body_is_fine_for_an_interpreter_that_does_not_care(self):
        """A `cat` of JSON, and SQL, arrive as data. Flagging them would be noise.

        `yq` is in that group too, which is not obvious: YAML is whitespace sensitive, but
        it fixes the root node's indentation from the first line, so a uniformly shifted
        document parses to the same value. Only `---` and `%YAML` need column zero.
        """
        for opener, body in (('cat > x.json <<FIN', '{"a": 1}'),
                             ('cat > x.sql <<FIN', 'select 1;'),
                             ('yq -o=json <<FIN', 'a: 1'),
                             ('R --vanilla <<FIN', 'print(1)')):
            with self.subTest(opener=opener):
                self.assertEqual([], broken(wdl([
                    'at_column_zero=1',
                    f'        {opener}',
                    f'        {body}',
                    'FIN',
                ])))

    def test_sensitivity_is_decided_by_the_command_not_the_delimiter(self):
        """`CODE`, `EOF` and `FIN` all name Python bodies somewhere in this repo."""
        self.assertEqual(1, len(broken(wdl([
            'at_column_zero=1',
            '        python3 - "$arg" > out.json <<FIN',
            '        import sys',
            'FIN',
        ]))))


class TestInlinePythonCompiles(unittest.TestCase):
    """Nothing else checks it: pyflakes cannot see Python embedded in a WDL string."""

    def test_a_syntax_error_is_reported(self):
        found = broken(wdl([
            "        python3 <<'CODE'",
            '        arguments = {',
            '            "action" action,',
            '        }',
            '        CODE',
        ]))
        self.assertEqual(1, len(found))
        self.assertIn('does not compile', found[0].message)

    def test_placeholders_do_not_look_like_syntax_errors(self):
        """Both interpolation styles, bare and quoted, as the WDLs here use them."""
        self.assertEqual([], checker.findings(wdl([
            "        python3 <<'CODE'",
            '        size = ~{superpartition_size}',
            '        mode = "~{mode}"',
            '        legacy = "${old_style}"',
            "        optional = \"~{default='' sample_map_path}\"",
            '        print(size, mode, legacy, optional)',
            '        CODE',
        ])))

    def test_an_indented_body_is_reported_once_not_twice(self):
        """IndentationError is a SyntaxError, so the compile check sees it too."""
        found = broken(wdl([
            'at_column_zero=1',
            "        python3 <<'CODE'",
            '        import sys',
            'CODE',
        ]))
        self.assertEqual(1, len(found), 'the indent and the IndentationError are one bug')
        self.assertIn('still indented', found[0].message)


class TestParsing(unittest.TestCase):

    def test_an_unterminated_command_block_is_still_checked(self):
        """womtool would reject the file, but silently checking nothing is worse."""
        text = '\n'.join(['task T {', '    command <<<', 'at_column_zero=1',
                          "        python3 <<'CODE'", '        import sys',
                          '        CODE'])
        self.assertEqual(1, len(broken(text)))

    def test_a_file_with_no_command_block_is_clean(self):
        self.assertEqual([], checker.findings('workflow W {\n    String x = "y"\n}\n'))

    def test_an_empty_command_block_is_clean(self):
        self.assertEqual([], checker.findings(wdl(['', '   '])))


class TestExitStatus(unittest.TestCase):
    """The script is a CI gate, so what it returns matters as much as what it prints."""

    @contextlib.contextmanager
    def written(self, block):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / 'T.wdl'
            path.write_text(wdl(block))
            yield str(path)

    def run_main(self, argv):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = checker.main(argv)
        return status, output.getvalue()

    def test_broken_fails(self):
        with self.written(['at_column_zero=1', "        python3 <<'CODE'",
                           '        import sys', '        CODE']) as path:
            status, output = self.run_main([path])
        self.assertEqual(1, status)
        self.assertIn('1 broken', output)

    def test_fragile_passes_by_default(self):
        with self.written(['echo hello']) as path:
            status, output = self.run_main([path])
        self.assertEqual(0, status)
        self.assertIn('0 broken, 1 fragile', output)

    def test_fragile_fails_under_strict(self):
        with self.written(['echo hello']) as path:
            status, _ = self.run_main(['--strict', path])
        self.assertEqual(1, status)

    def test_clean_passes(self):
        with self.written(['        echo hello']) as path:
            status, output = self.run_main([path])
        self.assertEqual(0, status)
        self.assertIn('0 broken, 0 fragile', output)


class TestEveryVariantstoreWdl(unittest.TestCase):
    """The point of the whole exercise: nothing in the tree we own is broken.

    Scoped to `scripts/variantstore` because that is what this team owns; the CNV and
    Mutect WDLs elsewhere in the repo have their own conventions and owners.
    """

    def test_no_variantstore_wdl_is_broken(self):
        root = checker.default_root()
        paths = checker.wdl_files([root])
        if not paths:
            # An explicit skip rather than a vacuous pass. A sweep that silently finds no
            # files is how a check stops guarding anything without ever failing, which is
            # the same shape as the bug this script exists to catch.
            self.skipTest(f'no WDL files under {root} in this test environment')
        broken = [f'{path}:{finding.line}: {finding.message}'
                  for path in paths
                  for finding in checker.findings(path.read_text())
                  if finding.severity == 'BROKEN']
        self.assertEqual([], broken, '\n'.join(broken))

    def test_the_sweep_actually_reaches_the_wdls(self):
        """Guards the guard: the assertion above is worthless if this is ever zero."""
        paths = checker.wdl_files([checker.default_root()])
        if not paths:
            self.skipTest('WDL tree not available in this test environment')
        self.assertGreater(len(paths), 20, 'the variantstore tree holds dozens of WDLs')


if __name__ == '__main__':
    unittest.main()
