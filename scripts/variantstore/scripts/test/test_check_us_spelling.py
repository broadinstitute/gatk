#!/usr/bin/env python3
"""Unit tests for the check-us-spelling helper.

check-us-spelling: disable -- these are fixtures, so they necessarily contain British forms.

The value of this tool is entirely in its precision. A spell checker that flags correct
prose gets switched off, so the false-positive cases below matter more than the true
positives: "analysis", "optimistic", "premise" and "exercise" are all correct US English
and a naive "-ise" rule reports every one of them.

The helper is deliberately extensionless -- see its docstring for why -- so it is loaded by
path rather than imported by name.
"""

import contextlib
import importlib.util
import io
import os
import pathlib
import tempfile
import unittest

_HELPER = pathlib.Path(__file__).resolve().parents[1] / 'check-us-spelling'
_spec = importlib.util.spec_from_loader(
    'check_us_spelling',
    importlib.machinery.SourceFileLoader('check_us_spelling', str(_HELPER)),
)
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)


class TestCorrectProseIsNotFlagged(unittest.TestCase):
    """Only the verbs shift: analyse -> analyze, but analysis stays put."""

    CORRECT = (
        'analysis', 'analyses', 'optimistic', 'premise', 'premises', 'exercise',
        'exercised', 'raise', 'raised', 'noise', 'promise', 'promised', 'precise',
        'concise', 'revised', 'advise', 'advised', 'surprise', 'supervise', 'comprise',
        'arise', 'arising', 'otherwise', 'likewise', 'wise', 'expertise', 'franchise',
        'merchandise', 'compromise', 'compromised', 'paradise', 'devise', 'improvise',
        'parameter', 'diameter', 'because', 'dense', 'phase', 'based', 'false', 'sense',
        'towards',
    )

    def test_none_are_flagged(self):
        for word in self.CORRECT:
            self.assertEqual([], list(checker.find_in_text(word)), msg=word)

    def test_a_paragraph_of_correct_prose_is_clean(self):
        text = ('The analysis was optimistic. Our premise: exercise raises noise, and we '
                'promise precise revised advice. Parameter and diameter arise otherwise.')
        self.assertEqual([], list(checker.find_in_text(text)))


class TestBritishFormsAreFlagged(unittest.TestCase):

    CASES = {
        'behaviour': 'behavior',
        'behavioural': 'behavioral',
        'recognised': 'recognized',
        'normalise': 'normalize',
        'normalisation': 'normalization',
        'prioritisation': 'prioritization',
        'analysed': 'analyzed',
        'analysing': 'analyzing',
        'colour': 'color',
        'centre': 'center',
        'defence': 'defense',
        'licence': 'license',
        'practised': 'practiced',
        'artefact': 'artifact',
        'judgement': 'judgment',
        'acknowledgement': 'acknowledgment',
        'fulfil': 'fulfill',
        'travelling': 'traveling',
        'modelling': 'modeling',
        'cancelled': 'canceled',
        'programme': 'program',
        'grey': 'gray',
        'whilst': 'while',
        'amongst': 'among',
        'learnt': 'learned',
        'sceptical': 'skeptical',
        'ageing': 'aging',
        'tokenise': 'tokenize',
        'criticised': 'criticized',
        'summarising': 'summarizing',
    }

    def test_each_is_flagged_with_the_right_suggestion(self):
        for uk, us in self.CASES.items():
            hits = list(checker.find_in_text(uk))
            self.assertEqual(1, len(hits), msg=uk)
            self.assertEqual((1, uk, us), hits[0], msg=uk)

    def test_the_words_that_slipped_into_this_repo(self):
        """Regression: these are the forms actually found during VS-1998."""
        for uk in ('behaviour', 'recognised'):
            self.assertTrue(list(checker.find_in_text(uk)), msg=uk)


class TestCasePreservation(unittest.TestCase):

    def test_lowercase(self):
        self.assertEqual('behavior', checker.fix_text('behaviour'))

    def test_capitalized(self):
        self.assertEqual('Behavior', checker.fix_text('Behaviour'))

    def test_all_caps(self):
        self.assertEqual('BEHAVIOR', checker.fix_text('BEHAVIOUR'))

    def test_mixed_case_degrades_to_lowercase_rather_than_guessing(self):
        self.assertEqual('behavior', checker.fix_text('bEhAvIoUr'))

    def test_surrounding_text_is_untouched(self):
        self.assertEqual('The behavior was recognized.',
                         checker.fix_text('The behaviour was recognised.'))


class TestWordBoundaries(unittest.TestCase):

    def test_substring_occurrences_are_not_matched(self):
        """`grey` must not fire inside a longer token, and vice versa."""
        for text in ('greyhound', 'centred_on_disk', 'xcentre', 'metres_per', 'ametre'):
            hits = [h for h in checker.find_in_text(text)]
            self.assertEqual([], hits, msg=text)

    def test_identifiers_with_underscores_are_out_of_scope(self):
        """Deliberate: --fix must never turn a spelling nit into broken code.

        Python's \\b treats underscore as a word character, so a snake_case identifier is
        not matched. Renaming one here but not at its definition in another file would
        break the build, so prose and comments are the target and identifiers are not.
        """
        self.assertEqual([], list(checker.find_in_text('n_behaviour_flags')))
        self.assertEqual([], list(checker.find_in_text('recognised_forms')))

    def test_a_bare_word_adjacent_to_punctuation_is_matched(self):
        for text in ('behaviour.', '(behaviour)', '"behaviour"', 'behaviour, colour'):
            self.assertTrue(list(checker.find_in_text(text)), msg=text)

    def test_longest_match_wins(self):
        hits = list(checker.find_in_text('normalisation'))
        self.assertEqual([(1, 'normalisation', 'normalization')], hits)


class TestLineNumbers(unittest.TestCase):

    def test_reported_per_line(self):
        text = 'clean line\nbehaviour here\nalso recognised\n'
        self.assertEqual(
            [(2, 'behaviour', 'behavior'), (3, 'recognised', 'recognized')],
            list(checker.find_in_text(text)))

    def test_multiple_hits_on_one_line(self):
        hits = list(checker.find_in_text('behaviour and colour'))
        self.assertEqual([(1, 'behaviour', 'behavior'), (1, 'colour', 'color')], hits)


class TestCommandLine(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_quietly(self, argv):
        """Invoke main() without its report reaching the test output."""
        with contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            return checker.main(argv)

    def write(self, name, content):
        path = os.path.join(self.tmp.name, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as handle:
            handle.write(content)
        return path

    def test_clean_file_exits_zero(self):
        path = self.write('clean.md', 'The analysis was optimistic.\n')
        self.assertEqual(0, self.run_quietly([path]))

    def test_dirty_file_exits_nonzero(self):
        """Non-zero so it works as a CI or pre-commit gate."""
        path = self.write('dirty.md', 'The behaviour was recognised.\n')
        self.assertEqual(1, self.run_quietly([path]))

    def test_fix_rewrites_and_exits_zero(self):
        path = self.write('dirty.md', 'The behaviour was recognised.\n')
        self.assertEqual(0, self.run_quietly(['--fix', path]))
        self.assertEqual('The behavior was recognized.\n', self.read(path))

    def read(self, path):
        with open(path) as handle:
            return handle.read()

    def test_fix_is_idempotent(self):
        path = self.write('dirty.md', 'Behaviour and colour.\n')
        self.run_quietly(['--fix', path])
        first = self.read(path)
        self.assertEqual(0, self.run_quietly(['--fix', path]))
        self.assertEqual(first, self.read(path))

    def test_directory_is_searched_recursively(self):
        self.write('sub/nested.md', 'behaviour\n')
        self.assertEqual(1, self.run_quietly([self.tmp.name]))

    def test_caches_and_binaries_are_skipped(self):
        self.write('__pycache__/x.md', 'behaviour\n')
        self.write('stale.pyc', 'behaviour\n')
        self.assertEqual(0, self.run_quietly([self.tmp.name]))

    def test_no_arguments_is_a_usage_error(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(2, checker.main([]))

    def test_unknown_flag_is_a_usage_error(self):
        path = self.write('clean.md', 'fine\n')
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(2, checker.main(['--nonsense', path]))

    def test_unreadable_file_is_skipped_not_fatal(self):
        path = os.path.join(self.tmp.name, 'binary.dat')
        with open(path, 'wb') as handle:
            handle.write(b'\xff\xfe\x00behaviour')
        self.assertEqual(0, self.run_quietly([path]))


class TestSuppression(unittest.TestCase):
    """Some files legitimately contain British forms and must be exemptable.

    Without this the tool cannot be put in CI: it would fail forever on its own word list
    and on this test file.
    """

    def test_file_marker_exempts_the_whole_file(self):
        text = 'x check-us-spelling: disable\nbehaviour and colour and recognised\n'
        self.assertEqual([], list(checker.find_in_text(text)))

    def test_line_marker_exempts_only_that_line(self):
        text = ('behaviour here\n'
                'colour here check-us-spelling: ok\n'
                'recognised here\n')
        hits = list(checker.find_in_text(text))
        self.assertEqual([(1, 'behaviour', 'behavior'), (3, 'recognised', 'recognized')],
                         hits)

    def test_fix_honors_the_file_marker(self):
        """An exempt file must never be rewritten, or the hatch only half works."""
        text = 'x check-us-spelling: disable\nbehaviour\n'
        self.assertEqual(text, checker.fix_text(text))

    def test_fix_honors_the_line_marker(self):
        text = 'behaviour\ncolour check-us-spelling: ok\n'
        self.assertEqual('behavior\ncolour check-us-spelling: ok\n',
                         checker.fix_text(text))

    def test_the_tool_and_its_test_are_exempt(self):
        """Both necessarily contain British forms; neither should ever be reported."""
        here = pathlib.Path(__file__).resolve()
        for path in (here, here.parents[1] / 'check-us-spelling'):
            self.assertEqual([], list(checker.find_in_text(path.read_text())),
                             msg=str(path))


class TestMappingHygiene(unittest.TestCase):

    def test_no_mapping_is_a_no_op(self):
        """A key equal to its value would report a word as needing itself."""
        for uk, us in checker.WORD_MAP.items():
            self.assertNotEqual(uk, us, msg=uk)

    def test_keys_are_lowercase(self):
        for uk in checker.WORD_MAP:
            self.assertEqual(uk, uk.lower(), msg=uk)

    def test_replacements_are_never_themselves_flagged(self):
        """Otherwise --fix would not converge."""
        for us in set(checker.WORD_MAP.values()):
            self.assertEqual([], list(checker.find_in_text(us)), msg=us)

    def test_deliberate_omissions_stay_omitted(self):
        """`analyses` is ambiguous and `towards` is style, not spelling.

        Both are documented as intentionally absent; re-adding either would produce false
        positives on correct prose.
        """
        for word in ('analyses', 'towards'):
            self.assertNotIn(word, checker.WORD_MAP, msg=word)


if __name__ == '__main__':
    unittest.main()
