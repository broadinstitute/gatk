#!/usr/bin/env python3
"""Unit tests for the jira-adf helper.

The helper is deliberately extensionless -- see its docstring for why -- so it is loaded by path
rather than imported by name.

The fixtures are synthetic. The guards matter more than the happy paths here: the helper exists to
make an edit that goes back to a shared JIRA ticket, so the tests concentrate on the cases where it
must refuse to write anything.
"""

import contextlib
import importlib.machinery
import importlib.util
import io
import json
import os
import pathlib
import tempfile
import unittest

_HELPER = pathlib.Path(__file__).resolve().parents[1] / 'jira-adf'
_spec = importlib.util.spec_from_loader(
    'jira_adf',
    importlib.machinery.SourceFileLoader('jira_adf', str(_HELPER)),
)
jira_adf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(jira_adf)


def text(s, marks=None):
    node = {'type': 'text', 'text': s}
    if marks:
        node['marks'] = marks
    return node


def para(*children, local_id=None):
    node = {'type': 'paragraph', 'content': list(children)}
    if local_id:
        node['attrs'] = {'localId': local_id}
    return node


def bullets(*items):
    return {'type': 'bulletList', 'attrs': {'localId': 'list0001'},
            'content': [{'type': 'listItem', 'attrs': {'localId': f'item{i:04d}'},
                         'content': [para(text(t))]} for i, t in enumerate(items)]}


def sample_doc():
    # The localId on block 1 deliberately contains "9419" while its visible text does not, so a
    # guard that matched against serialized JSON instead of visible text would wrongly pass.
    return {'type': 'doc', 'version': 1, 'content': [
        para(text('Intro paragraph.'), local_id='aaaa00000001'),
        para(text('Then two things to check:'), local_id='bbbb94190002'),
        bullets('PR #9419 has merged', 'VS-2021 has landed'),
        {'type': 'heading', 'attrs': {'level': 2, 'localId': 'cccc00000003'},
         'content': [text('Open items')]},
        bullets('No cluster can rotate zones', 'Rename the parameters', 'Pick the output name'),
        para(text('See '), {'type': 'inlineCard', 'attrs': {'url': 'https://example.invalid/browse/VS-1'}},
             text(' and '), text('code', [{'type': 'code'}]), text('.')),
        {'type': 'table', 'attrs': {'layout': 'default', 'localId': 'tbl00001'}, 'content': [
            {'type': 'tableRow', 'content': [
                {'type': 'tableCell', 'content': [para(text('Zone rotation'))]},
                {'type': 'tableCell', 'content': [para(text('VS-2021'))]},
            ]},
        ]},
    ]}


class HelperTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = pathlib.Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, obj):
        path = self.dir / name
        path.write_text(json.dumps(obj), encoding='utf-8')
        return str(path)

    def run_helper(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            jira_adf.main(list(argv))
        return out.getvalue()

    def assert_refuses(self, *argv, message):
        out_path = self.dir / 'out.json'
        if out_path.exists():
            out_path.unlink()
        with self.assertRaises(SystemExit) as cm:
            self.run_helper(*argv)
        self.assertIn(message, str(cm.exception))
        self.assertFalse(out_path.exists(), 'a refused edit must not write its output file')


class TestExtract(HelperTestCase):
    def test_extracts_description_from_saved_issue_result(self):
        doc = sample_doc()
        src = self.write('result.txt', {'key': 'VS-1', 'fields': {'description': doc}})
        out = str(self.dir / 'desc.json')
        self.run_helper('extract', src, '-o', out)
        self.assertEqual(json.loads(pathlib.Path(out).read_text()), doc)

    def test_accepts_a_bare_doc(self):
        doc = sample_doc()
        out = str(self.dir / 'desc.json')
        self.run_helper('extract', self.write('bare.json', doc), '-o', out)
        self.assertEqual(json.loads(pathlib.Path(out).read_text()), doc)

    def test_markdown_fetch_is_rejected_with_a_hint(self):
        src = self.write('md.txt', {'fields': {'description': '**TLDR** some markdown'}})
        with self.assertRaises(SystemExit) as cm:
            self.run_helper('extract', src, '-o', str(self.dir / 'x.json'))
        self.assertIn('responseContentFormat', str(cm.exception))


class TestOutline(HelperTestCase):
    def test_top_level_and_nested_paths(self):
        path = self.write('d.json', sample_doc())
        top = self.run_helper('outline', path).splitlines()
        self.assertEqual(len(top), 7)
        self.assertTrue(top[3].startswith('3 ') and 'heading2' in top[3] and 'Open items' in top[3])
        self.assertIn('https://example.invalid/browse/VS-1', top[5])
        nested = self.run_helper('outline', path, '--path', '4').splitlines()
        self.assertEqual([line.split()[0] for line in nested], ['4.0', '4.1', '4.2'])
        self.assertIn('No cluster', nested[0])


class TestDelete(HelperTestCase):
    def test_deletes_top_level_and_nested_blocks_by_original_numbering(self):
        path = self.write('d.json', sample_doc())
        out = str(self.dir / 'out.json')
        # Listed in ascending order on purpose: deleting 1 first would shift 2 onto the heading and
        # 4.0 onto the wrong list, so this only passes if deletions run in reverse order.
        self.run_helper('delete', path, '1:Then two', '2:#9419', '4.0:No cluster', '-o', out)
        result = json.loads(pathlib.Path(out).read_text())
        texts = [jira_adf.visible_text(b) for b in result['content']]
        self.assertEqual(len(texts), 5)
        self.assertEqual(texts[1], 'Open items')
        items = [jira_adf.visible_text(i) for i in result['content'][2]['content']]
        self.assertEqual(items, ['Rename the parameters', 'Pick the output name'])

    def test_guard_mismatch_refuses_and_writes_nothing(self):
        path = self.write('d.json', sample_doc())
        self.assert_refuses('delete', path, '0:Intro', '3:Something else', '-o', str(self.dir / 'out.json'),
                            message='does not contain')

    def test_guard_ignores_attributes(self):
        # Block 1's localId contains 9419 but its visible text does not.
        self.assertIn('9419', json.dumps(sample_doc()['content'][1]))
        path = self.write('d.json', sample_doc())
        self.assert_refuses('delete', path, '1:9419', '-o', str(self.dir / 'out.json'),
                            message='does not contain')

    def test_missing_path_refuses(self):
        path = self.write('d.json', sample_doc())
        self.assert_refuses('delete', path, '4.9', '-o', str(self.dir / 'out.json'),
                            message='does not exist')

    def test_nested_overlap_refuses(self):
        path = self.write('d.json', sample_doc())
        self.assert_refuses('delete', path, '4', '4.1', '-o', str(self.dir / 'out.json'),
                            message='is inside')


class TestReplace(HelperTestCase):
    def test_replaces_exactly_one_occurrence(self):
        path = self.write('d.json', sample_doc())
        out = str(self.dir / 'out.json')
        self.run_helper('replace', path, 'Zone rotation', 'Zone rotation and n2d defaults', '-o', out)
        result = json.loads(pathlib.Path(out).read_text())
        self.assertIn('Zone rotation and n2d defaults', jira_adf.visible_text(result['content'][6]))

    def test_ambiguous_text_refuses(self):
        # "VS-2021" is both a list item and a table cell.
        path = self.write('d.json', sample_doc())
        self.assert_refuses('replace', path, 'VS-2021', 'X', '-o', str(self.dir / 'out.json'),
                            message='occurs in 2 text nodes')

    def test_text_spanning_formatting_is_not_found(self):
        path = self.write('d.json', sample_doc())
        self.assert_refuses('replace', path, 'and code', 'X', '-o', str(self.dir / 'out.json'),
                            message='occurs in 0 text nodes')


class TestStripAndSplit(HelperTestCase):
    def test_strip_removes_every_local_id_and_nothing_else(self):
        doc = sample_doc()
        path = self.write('d.json', doc)
        out = self.dir / 'min.json'
        self.run_helper('strip', path, '-o', str(out))
        raw = out.read_text(encoding='utf-8')
        self.assertNotIn('localId', raw)
        self.assertNotIn('\n', raw)
        result = json.loads(raw)
        # Other attributes survive, and an attrs object left empty is dropped rather than kept as {}.
        self.assertEqual(result['content'][3]['attrs'], {'level': 2})
        self.assertEqual(result['content'][6]['attrs'], {'layout': 'default'})
        self.assertNotIn('attrs', result['content'][0])
        self.assertEqual(result['content'][5]['content'][1]['attrs']['url'],
                         'https://example.invalid/browse/VS-1')
        self.assertEqual(jira_adf.visible_text(result), jira_adf.visible_text(doc))

    def test_non_ascii_is_kept_literal(self):
        doc = {'type': 'doc', 'version': 1, 'content': [para(text('a — b'))]}
        out = self.dir / 'min.json'
        self.run_helper('strip', self.write('d.json', doc), '-o', str(out))
        self.assertIn('—', out.read_text(encoding='utf-8'))

    def test_split_parts_concatenate_to_the_input(self):
        doc = {'type': 'doc', 'version': 1,
               'content': [para(text(f'block {i} — café')) for i in range(200)]}
        src = self.dir / 'min.json'
        self.run_helper('strip', self.write('d.json', doc), '-o', str(src))
        prefix = str(self.dir / 'part')
        self.run_helper('split', str(src), '--max-chars', '997', '-o', prefix)
        parts = sorted(self.dir.glob('part.*'), key=lambda p: int(p.suffix[1:]))
        self.assertGreater(len(parts), 5)
        self.assertTrue(all(len(p.read_text(encoding='utf-8')) <= 997 for p in parts))
        joined = ''.join(p.read_text(encoding='utf-8') for p in parts)
        self.assertEqual(joined, src.read_text(encoding='utf-8'))


if __name__ == '__main__':
    unittest.main()
