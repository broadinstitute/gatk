"""
Unit tests for filter_set_site_stats_from_vds.py.

Covers the three behaviors called out by code review:
  1. Empty filter sets are rendered as "PASS".
  2. Non-empty filter sets are sorted alphabetically before being joined,
     so the grouping key is canonical regardless of set-iteration order.
  3. Sites that are monomorphic in the remaining samples (all-ref LGT) are
     excluded when require_non_ref=True.

Tests run with Hail in local mode and do not require a real VDS on disk or
a Dataproc cluster.  The entire suite is skipped gracefully when Hail is not
installed in the test environment.
"""

import unittest

try:
    import hail as hl
    HAIL_AVAILABLE = True
except ImportError:
    HAIL_AVAILABLE = False

if HAIL_AVAILABLE:
    from filter_set_site_stats_from_vds import filter_stats_from_rows_table


def _collect_as_dict(stats_table: 'hl.Table') -> dict:
    """Materialize a stats Table to a {filter_string: count} Python dict."""
    return {row.filters: row.n_sites for row in stats_table.collect()}


@unittest.skipUnless(HAIL_AVAILABLE, 'hail not installed – skipping')
class TestFilterStatsFromRowsTable(unittest.TestCase):
    """Tests for filter_stats_from_rows_table() using synthetic Hail Tables."""

    @classmethod
    def setUpClass(cls):
        hl.init(quiet=True, idempotent=True)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _make_rows_table(filter_sets: list) -> 'hl.Table':
        """Build a minimal Hail Table with a single ``filters: set<str>`` field."""
        return hl.Table.parallelize(
            [{'filters': fs} for fs in filter_sets],
            schema=hl.tstruct(filters=hl.tset(hl.tstr)),
        )

    # ------------------------------------------------------------------
    # filter_str rendering
    # ------------------------------------------------------------------

    def test_empty_filters_rendered_as_pass(self):
        """Sites with an empty filter set must appear under the key 'PASS'."""
        ht = self._make_rows_table([set(), set(), set()])
        result = _collect_as_dict(filter_stats_from_rows_table(ht))
        self.assertEqual(result, {'PASS': 3})

    def test_single_filter_preserved(self):
        """A single-element filter set is returned as the bare filter name."""
        ht = self._make_rows_table([{'LowQual'}, {'LowQual'}])
        result = _collect_as_dict(filter_stats_from_rows_table(ht))
        self.assertEqual(result, {'LowQual': 2})

    def test_multi_filter_sorted_alphabetically(self):
        """
        Multi-element filter sets are sorted alphabetically before joining.

        'EXCESS_ALLELES' < 'ExcessHet' in Unicode order (uppercase 'X' < lowercase 'x'),
        so the canonical key is 'EXCESS_ALLELES, ExcessHet' regardless of
        the iteration order of the input Python set.
        """
        # Pass set in non-deterministic order; sorted output must be stable.
        ht = self._make_rows_table([{'ExcessHet', 'EXCESS_ALLELES'}])
        result = _collect_as_dict(filter_stats_from_rows_table(ht))
        self.assertEqual(list(result.keys()), ['EXCESS_ALLELES, ExcessHet'])
        self.assertEqual(result['EXCESS_ALLELES, ExcessHet'], 1)

    def test_multi_filter_sort_with_lowercase_first(self):
        """
        'LowQual' < 'NO_HQ_GENOTYPES' alphabetically, so sorted key is
        'LowQual, NO_HQ_GENOTYPES' (matches the original BigQuery-derived counts).
        """
        ht = self._make_rows_table([{'NO_HQ_GENOTYPES', 'LowQual'}])
        result = _collect_as_dict(filter_stats_from_rows_table(ht))
        self.assertEqual(list(result.keys()), ['LowQual, NO_HQ_GENOTYPES'])

    # ------------------------------------------------------------------
    # Grouping and counting
    # ------------------------------------------------------------------

    def test_counts_grouped_correctly(self):
        """Sites with the same canonical filter string are counted together."""
        ht = self._make_rows_table([
            set(),                              # → PASS
            {'LowQual'},                        # → LowQual
            {'LowQual'},                        # → LowQual
            {'LowQual'},                        # → LowQual  (count=3)
            {'NO_HQ_GENOTYPES', 'LowQual'},     # → LowQual, NO_HQ_GENOTYPES
            {'LowQual', 'NO_HQ_GENOTYPES'},     # → same key after sort (count=2)
        ])
        result = _collect_as_dict(filter_stats_from_rows_table(ht))
        self.assertEqual(result.get('PASS'), 1)
        self.assertEqual(result.get('LowQual'), 3)
        self.assertEqual(result.get('LowQual, NO_HQ_GENOTYPES'), 2)
        # No other keys should appear.
        self.assertEqual(set(result.keys()), {'PASS', 'LowQual', 'LowQual, NO_HQ_GENOTYPES'})

    def test_mixed_pass_and_filtered(self):
        """PASS sites and filtered sites coexist correctly in the output."""
        ht = self._make_rows_table([set(), {'ExcessHet'}])
        result = _collect_as_dict(filter_stats_from_rows_table(ht))
        self.assertEqual(result.get('PASS'), 1)
        self.assertEqual(result.get('ExcessHet'), 1)


@unittest.skipUnless(HAIL_AVAILABLE, 'hail not installed – skipping')
class TestMonomorphicRowFiltering(unittest.TestCase):
    """
    Tests for the require_non_ref guard that drops rows where no sample
    carries a non-reference genotype.

    These tests apply the filter expression directly against a synthetic
    MatrixTable rather than calling compute_filter_set_site_stats(), because
    that function requires a real VDS on disk.  The filter logic inside
    compute_filter_set_site_stats() is a single hl.agg.any(mt.LGT.is_non_ref())
    call — the same expression exercised here — so the coverage gap is minimal.
    """

    @classmethod
    def setUpClass(cls):
        hl.init(quiet=True, idempotent=True)

    def _make_variant_mt(self) -> 'hl.MatrixTable':
        """
        Build a 3-row × 2-column MatrixTable with:
          row 0: sample 0 → 0/1 (non-ref), sample 1 → 0/0 (ref)  → KEPT
          row 1: sample 0 → 0/0,           sample 1 → 0/0         → DROPPED
          row 2: sample 0 → 0/0,           sample 1 → 0/1         → KEPT

        All rows have an empty filters set (PASS).
        """
        mt = hl.utils.range_matrix_table(n_rows=3, n_cols=2)
        mt = mt.annotate_rows(filters=hl.empty_set(hl.tstr))
        mt = mt.annotate_entries(
            LGT=hl.case()
                .when((mt.row_idx == 0) & (mt.col_idx == 0), hl.call(0, 1))
                .when((mt.row_idx == 2) & (mt.col_idx == 1), hl.call(0, 1))
                .default(hl.call(0, 0))
        )
        return mt

    def test_monomorphic_row_is_excluded(self):
        """
        Row 1 (all-ref) must be dropped; only the 2 non-ref rows should be
        counted under 'PASS'.
        """
        mt = self._make_variant_mt()
        filtered_mt = mt.filter_rows(hl.agg.any(mt.LGT.is_non_ref()))
        result = _collect_as_dict(filter_stats_from_rows_table(filtered_mt.rows()))
        self.assertEqual(result.get('PASS'), 2,
                         "Expected exactly 2 PASS sites after dropping the all-ref row")

    def test_without_filter_all_rows_present(self):
        """
        Without the non-ref filter, all 3 rows (including the monomorphic one)
        should be counted.
        """
        mt = self._make_variant_mt()
        result = _collect_as_dict(filter_stats_from_rows_table(mt.rows()))
        self.assertEqual(result.get('PASS'), 3,
                         "Expected all 3 rows when monomorphic filter is NOT applied")


if __name__ == '__main__':
    unittest.main()

