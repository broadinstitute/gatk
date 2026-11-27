import re
import unittest
from pathlib import Path

import populate_alt_allele_table as paat

TEMPLATE = (Path(__file__).resolve().parent.parent / 'alt_allele_positions.sql').read_text()


def _branches(sql):
    """The SELECT branches of the UNION ALL, one per alt_allele row shape."""
    return [b for b in re.split(r'\bunion all\b', sql, flags=re.IGNORECASE) if b.strip()]


class TestAltAllelePositionsSql(unittest.TestCase):

    def test_without_vrs_selects_no_id_column(self):
        # A dataset without VRS IDs has no vet.vrs_allele_ids and no alt_allele.vrs_allele_id, so the query must not
        # name either, or it fails on the vet table or inserts one column too many.
        sql = paat.alt_allele_positions_sql(TEMPLATE, include_vrs_allele_id=False)
        self.assertNotIn('vrs', sql)
        self.assertNotIn('$', sql)

    def test_with_vrs_reads_the_id_at_the_rows_alt_offset(self):
        # vet.vrs_allele_ids is aligned with SPLIT(alt, ','), so each branch must read the ID at the offset it reads
        # the ALT from, which is allele_pos - 1. A 1/2 call's second row is where an off-by-one would land.
        sql = paat.alt_allele_positions_sql(TEMPLATE, include_vrs_allele_id=True)
        self.assertNotIn('$', sql)
        branches = _branches(sql)
        self.assertEqual(len(branches), 3)
        seen = []
        for branch in branches:
            allele_pos = int(re.search(r'(\d+) as allele_pos', branch).group(1))
            alt_offsets = {int(n) for n in re.findall(r"SPLIT\(alt,','\)\[OFFSET\((\d+)\)\]", branch)}
            id_offsets = [int(n) for n in re.findall(r"NULLIF\(vrs_allele_ids\[OFFSET\((\d+)\)\], ''\)", branch)]
            self.assertEqual(alt_offsets, {allele_pos - 1}, branch)
            self.assertEqual(id_offsets, [allele_pos - 1], branch)
            seen.append(allele_pos)
        self.assertEqual(seen, [1, 1, 2])

    def test_with_vrs_the_id_is_the_last_column(self):
        # The INSERT into alt_allele is positional, and vrs_allele_id is the table's last column.
        sql = paat.alt_allele_positions_sql(TEMPLATE, include_vrs_allele_id=True)
        for branch in _branches(sql):
            self.assertRegex(branch, r"as vrs_allele_id\s+from position[12]")


if __name__ == '__main__':
    unittest.main()
