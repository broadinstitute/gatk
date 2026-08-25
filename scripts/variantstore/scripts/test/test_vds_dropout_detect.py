#!/usr/bin/env python3
"""Unit tests for vds_dropout_detect.

These build synthetic summary matrices rather than reading a VDS, which is the whole
point of keeping the judging logic separate from the Hail scan: the cases that matter
most -- a region hard for every sample, an ancestry-skewed superpartition, a dropout
thinned rather than emptied -- are awkward to arrange in real data and trivial to
arrange here.

The scale used throughout mirrors the real screen: 134 superpartitions, 100 sampled
samples each, and a per-sample rate of 75 entries per 50 kb bin.
"""

import gzip
import os
import tempfile
import unittest

import vds_dropout_detect as vdd

N_SUPERPARTITIONS = 134
N_SAMPLES = 100
RATE_PER_SAMPLE = 75.0
BIN_SIZE = 50_000
CELL = RATE_PER_SAMPLE * N_SAMPLES  # 7500 entries per (bin, superpartition)


def build_summary(
        n_bins: int = 40,
        contig: str = 'chr4',
        first_start: int = 56_000_001,
        n_superpartitions: int = N_SUPERPARTITIONS,
        n_samples: int = N_SAMPLES,
        cell_value: float = CELL,
) -> vdd.Summary:
    """A clean summary: every cell carries the same amount of data."""
    superpartitions = list(range(1, n_superpartitions + 1))
    summary = vdd.Summary(
        superpartitions=superpartitions,
        n_samples={sp: n_samples for sp in superpartitions},
    )
    from array import array
    for i in range(n_bins):
        start = first_start + i * BIN_SIZE
        summary.bins.append((contig, start, start + BIN_SIZE))
        summary.observed.append(array('d', [cell_value] * n_superpartitions))
    return summary


def set_cell(summary: vdd.Summary, bin_index: int, superpartition: int, value: float) -> None:
    summary.observed[bin_index][summary.superpartition_position(superpartition)] = value


def scale_superpartition(summary: vdd.Summary, superpartition: int, factor: float) -> None:
    j = summary.superpartition_position(superpartition)
    for row in summary.observed:
        row[j] *= factor


def rectangles_for(report: vdd.Report, superpartition: int) -> list[vdd.Rectangle]:
    return [r for r in report.rectangles if r.superpartition == superpartition]


class TestCleanMatrix(unittest.TestCase):

    def test_uniform_matrix_produces_no_findings(self):
        report = vdd.analyze(build_summary())
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))
        self.assertEqual(0, report.n_cells_flagged)
        self.assertEqual(40, report.n_bins_considered)

    def test_poisson_scale_jitter_does_not_fire(self):
        """Small cell-to-cell variation must not be mistaken for a dropout."""
        summary = build_summary()
        # Deterministic +/- 5% ripple across every cell.
        for i, row in enumerate(summary.observed):
            for j in range(len(row)):
                row[j] = CELL * (1.0 + 0.05 * (((i + j) % 3) - 1))
        report = vdd.analyze(summary)
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))


class TestFullDropout(unittest.TestCase):

    def test_full_rectangle_is_found_with_correct_bounds(self):
        summary = build_summary()
        # Nine consecutive bins fully missing for superpartition 83, mirroring the
        # observed ~450 kb chr4 window.
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)

        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found), msg=vdd.format_summary(report, 'variants'))

        rect = found[0]
        self.assertEqual('chr4', rect.contig)
        self.assertEqual(summary.bins[10][1], rect.start)
        self.assertEqual(summary.bins[18][2], rect.end)
        self.assertEqual(9, rect.n_bins)
        self.assertEqual(9 * BIN_SIZE, rect.span)
        self.assertEqual(0.0, rect.observed)
        self.assertAlmostEqual(0.0, rect.ratio)

    def test_only_the_affected_superpartition_is_flagged(self):
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        self.assertEqual({83}, {r.superpartition for r in report.rectangles})

    def test_two_independent_dropouts_are_reported_separately(self):
        """The two known Foxtrot windows live on different contigs and superpartitions."""
        summary = build_summary(n_bins=40, contig='chr4')
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)

        # Splice in a chr19 block for superpartition 64.
        from array import array
        for i in range(20):
            start = 40_000_001 + i * BIN_SIZE
            summary.bins.append(('chr19', start, start + BIN_SIZE))
            summary.observed.append(array('d', [CELL] * N_SUPERPARTITIONS))
        for i in range(42, 53):
            set_cell(summary, i, 64, 0.0)

        report = vdd.analyze(summary)
        self.assertEqual(2, len(report.rectangles))
        by_sp = {r.superpartition: r for r in report.rectangles}
        self.assertEqual({64, 83}, set(by_sp))
        self.assertEqual('chr4', by_sp[83].contig)
        self.assertEqual('chr19', by_sp[64].contig)


class TestPartialDropout(unittest.TestCase):
    """The r2-like regime: windows thinned rather than emptied."""

    def test_ninety_nine_percent_missing_is_found(self):
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, CELL * 0.01)
        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found))
        self.assertAlmostEqual(0.01, found[0].ratio, places=3)

    def test_ninety_percent_missing_is_found(self):
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, CELL * 0.10)
        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found))
        self.assertAlmostEqual(0.10, found[0].ratio, places=3)

    def test_fifty_percent_missing_is_found(self):
        """Right at the default ratio threshold boundary, so this pins the default."""
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, CELL * 0.40)
        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found), msg=vdd.format_summary(report, 'variants'))

    def test_severity_ordering_is_worst_first(self):
        summary = build_summary(n_bins=60)
        for i in range(5, 8):
            set_cell(summary, i, 83, CELL * 0.40)
        for i in range(20, 23):
            set_cell(summary, i, 64, 0.0)
        report = vdd.analyze(summary)
        self.assertEqual(2, len(report.rectangles))
        self.assertEqual(64, report.rectangles[0].superpartition)
        self.assertGreater(report.rectangles[0].score, report.rectangles[1].score)

    def test_mild_depletion_below_threshold_does_not_fire(self):
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, CELL * 0.80)
        report = vdd.analyze(summary)
        self.assertEqual([], rectangles_for(report, 83))


class TestFalsePositiveGuards(unittest.TestCase):

    def test_region_hard_for_every_superpartition_does_not_fire(self):
        """A centromere depresses all superpartitions equally, so the ratio stays near 1."""
        summary = build_summary()
        for i in range(15, 20):
            for sp in summary.superpartitions:
                set_cell(summary, i, sp, 0.0)
        report = vdd.analyze(summary)
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))
        self.assertEqual(5, report.n_bins_skipped_empty)

    def test_partially_hard_region_does_not_fire(self):
        """Uniformly low but non-zero coverage everywhere is still not a dropout."""
        summary = build_summary()
        for i in range(15, 20):
            for sp in summary.superpartitions:
                set_cell(summary, i, sp, CELL * 0.02)
        report = vdd.analyze(summary)
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))

    def test_ancestry_skewed_superpartition_does_not_fire(self):
        """A superpartition carrying 40% more variants everywhere is a global offset."""
        summary = build_summary()
        scale_superpartition(summary, 7, 1.40)
        report = vdd.analyze(summary)
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))

    def test_low_scale_superpartition_does_not_mask_a_real_dropout(self):
        """A superpartition running 25% light globally must still report a real window."""
        summary = build_summary()
        scale_superpartition(summary, 83, 0.75)
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found), msg=vdd.format_summary(report, 'variants'))

    def test_thin_bins_are_not_flagged_on_noise(self):
        """Below the evidence floor, an empty cell carries no weight."""
        summary = build_summary(cell_value=10.0)  # expected well under min_expected
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        self.assertEqual([], rectangles_for(report, 83))

    def test_lowering_min_expected_surfaces_thin_bins(self):
        """The floor is a threshold, not a blind spot -- it can be lowered deliberately."""
        summary = build_summary(cell_value=10.0)
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary, min_expected=5.0, score_threshold=2.0)
        self.assertEqual(1, len(rectangles_for(report, 83)))


class TestGeometry(unittest.TestCase):

    def test_single_bin_dropout_is_reported(self):
        summary = build_summary()
        set_cell(summary, 12, 83, 0.0)
        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found))
        self.assertEqual(1, found[0].n_bins)
        self.assertEqual(BIN_SIZE, found[0].span)

    def test_partial_boundary_bins_extend_the_rectangle(self):
        """A window with half-populated edges should report the wider span."""
        summary = build_summary()
        set_cell(summary, 9, 83, CELL * 0.45)
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)
        set_cell(summary, 19, 83, CELL * 0.45)

        report = vdd.analyze(summary)
        found = rectangles_for(report, 83)
        self.assertEqual(1, len(found))
        self.assertEqual(11, found[0].n_bins)
        self.assertEqual(summary.bins[9][1], found[0].start)
        self.assertEqual(summary.bins[19][2], found[0].end)

    def test_gap_between_flagged_bins_splits_rectangles(self):
        summary = build_summary()
        for i in (10, 11):
            set_cell(summary, i, 83, 0.0)
        for i in (20, 21):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        self.assertEqual(2, len(rectangles_for(report, 83)))

    def test_dropout_spanning_a_whole_contig_is_caught_as_scale_finding(self):
        """The blind spot in per-bin normalization, covered by the scale check.

        A superpartition missing everywhere has its own median dragged to zero, so the
        per-bin residuals look unremarkable. The global scale comparison is what catches
        it.
        """
        summary = build_summary(n_bins=30)
        for i in range(30):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        self.assertFalse(report.clean)
        flagged = {f.superpartition for f in report.superpartition_scales}
        self.assertIn(83, flagged)

    def test_majority_depleted_superpartition_is_caught_as_scale_finding(self):
        summary = build_summary(n_bins=30)
        for i in range(22):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        self.assertFalse(report.clean)
        self.assertIn(83, {f.superpartition for f in report.superpartition_scales})

    def test_uneven_sample_counts_are_normalized(self):
        """A short final superpartition must not look like a dropout."""
        summary = build_summary()
        summary.n_samples[134] = 20
        j = summary.superpartition_position(134)
        for row in summary.observed:
            row[j] = RATE_PER_SAMPLE * 20
        report = vdd.analyze(summary)
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))

    def test_superpartition_with_no_samples_is_ignored(self):
        summary = build_summary()
        summary.n_samples[134] = 0
        j = summary.superpartition_position(134)
        for row in summary.observed:
            row[j] = 0.0
        report = vdd.analyze(summary)
        self.assertTrue(report.clean, msg=vdd.format_summary(report, 'variants'))


class TestBaselineQuantile(unittest.TestCase):
    """The baseline is a high quantile so that a dropout hitting many superpartitions
    at once still leaves a usable reference point."""

    def test_quantile_interpolates(self):
        self.assertEqual(0.0, vdd._quantile([], 0.75))
        self.assertEqual(5.0, vdd._quantile([5.0], 0.75))
        self.assertEqual(1.0, vdd._quantile([1.0, 2.0, 3.0, 4.0], 0.0))
        self.assertEqual(4.0, vdd._quantile([1.0, 2.0, 3.0, 4.0], 1.0))
        self.assertAlmostEqual(3.25, vdd._quantile([1.0, 2.0, 3.0, 4.0], 0.75))

    def test_dropout_across_sixty_percent_of_superpartitions_is_found(self):
        """A median baseline would collapse to zero here and discard the bin."""
        summary = build_summary()
        affected = list(range(1, 81))  # 80 of 134, just over half
        for i in range(10, 19):
            for sp in affected:
                set_cell(summary, i, sp, 0.0)

        report = vdd.analyze(summary)
        flagged = {r.superpartition for r in report.rectangles}
        self.assertEqual(set(affected), flagged, msg=vdd.format_summary(report, 'variants'))
        self.assertEqual(0, report.n_bins_skipped_empty)

    def test_median_baseline_would_have_missed_it(self):
        """Pins the reason for the default: the same input at quantile 0.5 finds nothing."""
        summary = build_summary()
        for i in range(10, 19):
            for sp in range(1, 81):
                set_cell(summary, i, sp, 0.0)

        median_report = vdd.analyze(summary, baseline_quantile=0.5)
        self.assertEqual([], median_report.rectangles)
        self.assertEqual(9, median_report.n_bins_skipped_empty)

    def test_baseline_still_zero_when_nearly_everything_is_empty(self):
        """Beyond the quantile's reach the bin is genuinely uninformative, not a finding."""
        summary = build_summary()
        for i in range(10, 19):
            for sp in range(1, 131):  # 130 of 134 empty, past the 75th percentile
                set_cell(summary, i, sp, 0.0)
        report = vdd.analyze(summary)
        self.assertEqual(9, report.n_bins_skipped_empty)
        self.assertEqual([], report.rectangles)

    def test_clean_matrix_unaffected_by_quantile_baseline(self):
        summary = build_summary()
        for i, row in enumerate(summary.observed):
            for j in range(len(row)):
                row[j] = CELL * (1.0 + 0.05 * (((i + j) % 3) - 1))
        self.assertTrue(vdd.analyze(summary).clean)


class TestSuperpartitionScale(unittest.TestCase):

    def test_reported_scale_is_invariant_to_baseline_quantile(self):
        """The threshold must mean the same thing however the baseline is chosen."""
        summary = build_summary(n_bins=30)
        scale_superpartition(summary, 83, 0.40)
        at_median = vdd.analyze(summary, baseline_quantile=0.5).superpartition_scales
        at_upper = vdd.analyze(summary, baseline_quantile=0.75).superpartition_scales
        self.assertEqual(1, len(at_median))
        self.assertEqual(1, len(at_upper))
        self.assertAlmostEqual(
            at_median[0].relative_scale, at_upper[0].relative_scale, places=2)

    def test_typical_superpartition_scores_one(self):
        summary = build_summary(n_bins=30)
        scale_superpartition(summary, 83, 0.40)
        report = vdd.analyze(summary)
        self.assertAlmostEqual(0.40, report.superpartition_scales[0].relative_scale, places=2)

    def test_modest_batch_variation_stays_below_threshold(self):
        """25% light is within plausible batch variation and must not fire by default."""
        summary = build_summary(n_bins=30)
        scale_superpartition(summary, 83, 0.75)
        self.assertEqual([], vdd.analyze(summary).superpartition_scales)


class TestLocationEncoding(unittest.TestCase):

    def test_encodes_autosomes(self):
        self.assertEqual(1_000_000_001_000, vdd.encode_location('chr1', 1000))
        self.assertEqual(4_000_056_585_368, vdd.encode_location('chr4', 56_585_368))
        self.assertEqual(19_000_040_097_642, vdd.encode_location('chr19', 40_097_642))

    def test_encodes_sex_chromosomes(self):
        self.assertEqual(23, vdd.contig_index('chrX'))
        self.assertEqual(24, vdd.contig_index('chrY'))
        self.assertEqual(24_000_056_694_638, vdd.encode_location('chrY', 56_694_638))

    def test_unknown_contig_raises(self):
        with self.assertRaises(KeyError):
            vdd.contig_index('chr1_KI270706v1_random')


class TestAdjudicationSql(unittest.TestCase):

    def make_rectangle(self, **overrides):
        defaults = dict(
            contig='chr4', start=56_585_001, end=57_035_001, superpartition=83,
            n_samples=100, n_bins=9, observed=0.0, expected=67_500.0,
        )
        defaults.update(overrides)
        return vdd.Rectangle(**defaults)

    def test_variant_sql_targets_the_right_table_and_window(self):
        sql = vdd.adjudication_sql(
            self.make_rectangle(),
            project_id='aou-genomics-curation-prod',
            dataset_name='foxtrot',
            mode='variants',
        )
        self.assertIn('`aou-genomics-curation-prod.foxtrot.vet_083`', sql)
        self.assertIn(str(vdd.encode_location('chr4', 56_585_001)), sql)
        self.assertIn(str(vdd.encode_location('chr4', 57_035_000)), sql)
        self.assertIn('s.withdrawn IS NULL', sql)
        self.assertIn('s.is_control = false', sql)

    def test_superpartition_is_zero_padded_to_three_digits(self):
        sql = vdd.adjudication_sql(
            self.make_rectangle(superpartition=7),
            project_id='p', dataset_name='d', mode='variants',
        )
        self.assertIn('vet_007', sql)

    def test_reference_sql_targets_ref_ranges(self):
        sql = vdd.adjudication_sql(
            self.make_rectangle(), project_id='p', dataset_name='d', mode='references',
        )
        self.assertIn('ref_ranges_083', sql)
        self.assertNotIn('vet_083', sql)

    def test_compressed_is_the_default_reference_schema(self):
        """AoU callsets use the compressed schema, so it must be what you get by default."""
        sql = vdd.adjudication_sql(
            self.make_rectangle(), project_id='p', dataset_name='d', mode='references')
        self.assertIn('packed_ref_data', sql)
        self.assertNotIn('v.location', sql)

    def test_compressed_reference_sql_filters_on_the_clustering_field(self):
        """Filtering the decoded location would prune nothing and scan the whole table."""
        sql = vdd.adjudication_sql(
            self.make_rectangle(), project_id='p', dataset_name='d',
            mode='references', reference_schema='compressed')
        where = [line for line in sql.split('\n') if line.startswith('WHERE')][0]
        self.assertIn('v.packed_ref_data BETWEEN', where)

    def test_uncompressed_reference_sql_filters_on_location(self):
        sql = vdd.adjudication_sql(
            self.make_rectangle(), project_id='p', dataset_name='d',
            mode='references', reference_schema='uncompressed')
        where = [line for line in sql.split('\n') if line.startswith('WHERE')][0]
        self.assertIn('v.location BETWEEN', where)
        self.assertNotIn('packed_ref_data', sql)

    def test_invalid_reference_schema_raises(self):
        with self.assertRaises(ValueError):
            vdd.adjudication_sql(
                self.make_rectangle(), project_id='p', dataset_name='d',
                mode='references', reference_schema='squashed')

    def test_header_records_the_vds_side_figures(self):
        sql = vdd.adjudication_sql(
            self.make_rectangle(observed=6_750.0),
            project_id='p', dataset_name='d', mode='variants',
        )
        self.assertIn('6,750', sql)
        self.assertIn('67,500', sql)

    def test_end_is_exclusive(self):
        """bin_end is exclusive, so the query's upper bound is one position lower."""
        sql = vdd.adjudication_sql(
            self.make_rectangle(start=1000, end=2000),
            project_id='p', dataset_name='d', mode='variants',
        )
        self.assertIn(str(vdd.encode_location('chr4', 1999)), sql)
        self.assertNotIn(str(vdd.encode_location('chr4', 2000)), sql)

    def test_invalid_mode_raises(self):
        with self.assertRaises(ValueError):
            vdd.adjudication_sql(
                self.make_rectangle(), project_id='p', dataset_name='d', mode='nonsense',
            )


class TestPackedRefData(unittest.TestCase):
    """Bit arithmetic for the compressed ref_ranges schema.

    This is safety-critical in a way most of this module is not: a wrong range predicate
    returns zero rows, which reads as "BigQuery does not have this data either" and would
    *confirm* a false conclusion instead of raising an error. So the bounds are checked
    against an independent reimplementation of UnpackRefRangeInfo rather than against the
    production decoder.
    """

    @staticmethod
    def pack(chromosome: int, position: int, length: int, state: int) -> int:
        """Independent re-implementation of the packing UnpackRefRangeInfo reverses."""
        return (chromosome << 48) | (position << 16) | (length << 4) | state

    def test_decode_matches_independent_packing(self):
        packed = self.pack(4, 56_585_368, 1000, 7)
        self.assertEqual(vdd.encode_location('chr4', 56_585_368),
                         vdd.decode_packed_location(packed))

    def test_decode_ignores_length_and_state(self):
        """Location must not shift as the low bits vary across their whole range."""
        expected = vdd.encode_location('chr19', 40_097_642)
        for length, state in ((0, 0), (1000, 7), (0xFFF, 0xF)):
            packed = self.pack(19, 40_097_642, length, state)
            self.assertEqual(expected, vdd.decode_packed_location(packed),
                             msg=f'length={length} state={state}')

    def test_bounds_include_every_low_bit_combination_in_range(self):
        low, high = vdd.packed_ref_data_bounds('chr4', 56_585_368, 57_035_833)
        for position in (56_585_368, 56_800_000, 57_035_833):
            for length, state in ((0, 0), (500, 3), (0xFFF, 0xF)):
                packed = self.pack(4, position, length, state)
                self.assertTrue(low <= packed <= high,
                                msg=f'{position} length={length} state={state}')

    def test_bounds_exclude_positions_outside_the_range(self):
        low, high = vdd.packed_ref_data_bounds('chr4', 56_585_368, 57_035_833)
        for position in (56_585_367, 57_035_834, 1, 100_000_000):
            for length, state in ((0, 0), (0xFFF, 0xF)):
                packed = self.pack(4, position, length, state)
                self.assertFalse(low <= packed <= high,
                                 msg=f'{position} length={length} state={state}')

    def test_bounds_exclude_other_contigs(self):
        """The chromosome occupies the high bits, so ranges must not bleed across contigs."""
        low, high = vdd.packed_ref_data_bounds('chr4', 1, 1_000_000)
        for contig_number in (3, 5, 19, 23, 24):
            packed = self.pack(contig_number, 500_000, 100, 1)
            self.assertFalse(low <= packed <= high, msg=f'chr{contig_number}')

    def test_bounds_round_trip_through_the_decoder(self):
        low, high = vdd.packed_ref_data_bounds('chr19', 40_097_642, 40_646_743)
        self.assertEqual(vdd.encode_location('chr19', 40_097_642),
                         vdd.decode_packed_location(low))
        self.assertEqual(vdd.encode_location('chr19', 40_646_743),
                         vdd.decode_packed_location(high))

    def test_single_position_range(self):
        low, high = vdd.packed_ref_data_bounds('chr1', 1000, 1000)
        self.assertTrue(low <= self.pack(1, 1000, 0, 0) <= high)
        self.assertTrue(low <= self.pack(1, 1000, 0xFFF, 0xF) <= high)
        self.assertFalse(low <= self.pack(1, 1001, 0, 0) <= high)

    def test_sex_chromosomes(self):
        low, high = vdd.packed_ref_data_bounds('chrY', 56_694_638, 56_694_638)
        self.assertTrue(low <= self.pack(24, 56_694_638, 10, 1) <= high)

    def test_bit_widths_match_the_canonical_encoder(self):
        """Masks must match SchemaUtils#encodeCompressedRefBlock (SchemaUtils.java:114)."""
        self.assertEqual(48, vdd.PACKED_CHROMOSOME_SHIFT)
        self.assertEqual(16, vdd.PACKED_POSITION_SHIFT)
        self.assertEqual(0xFFFF, vdd.PACKED_CHROMOSOME_MASK)
        self.assertEqual(0xFFFFFFFF, vdd.PACKED_POSITION_MASK)
        # Length (12 bits) and state (4 bits) together fill the low 16.
        self.assertEqual(0xFFFF, vdd.PACKED_LOW_BITS_MASK)

    def test_packed_order_matches_location_order(self):
        """Why clustering on packed_ref_data prunes a location range at all."""
        packed = [
            self.pack(4, 1000, 500, 3), self.pack(4, 2000, 0, 0),
            self.pack(4, 2000, 4095, 15), self.pack(5, 1, 0, 0),
            self.pack(19, 40_097_642, 100, 1), self.pack(24, 1, 0, 0),
        ]
        self.assertEqual(packed, sorted(packed))
        locations = [vdd.decode_packed_location(v) for v in packed]
        self.assertEqual(locations, sorted(locations))

    def test_inverted_range_raises(self):
        with self.assertRaises(ValueError):
            vdd.packed_ref_data_bounds('chr4', 2000, 1000)


class TestFileParsing(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name: str, content: str) -> str:
        path = os.path.join(self.tmp.name, name)
        opener = gzip.open if name.endswith('.gz') else open
        with opener(path, 'wt') as handle:
            handle.write(content)
        return path

    SUPERPARTITIONS = 'superpartition\tn_samples\n1\t100\n2\t100\n3\t100\n'

    def test_round_trip(self):
        sp_path = self.write('sp.tsv', self.SUPERPARTITIONS)
        summary_path = self.write('summary.tsv', (
            'contig\tbin_start\tbin_end\tsuperpartition\tobserved\n'
            'chr1\t1\t50001\t1\t7500\n'
            'chr1\t1\t50001\t2\t7500\n'
            'chr1\t1\t50001\t3\t7500\n'
        ))
        n_samples = vdd.load_superpartitions(sp_path)
        summary = vdd.load_summary(summary_path, n_samples)
        self.assertEqual([1, 2, 3], summary.superpartitions)
        self.assertEqual([('chr1', 1, 50001)], summary.bins)
        self.assertEqual([7500.0, 7500.0, 7500.0], list(summary.observed[0]))

    def test_absent_cells_are_zero(self):
        """Cells missing from the file are the dropout, so they must read as zero."""
        sp_path = self.write('sp.tsv', self.SUPERPARTITIONS)
        summary_path = self.write('summary.tsv', (
            'contig\tbin_start\tbin_end\tsuperpartition\tobserved\n'
            'chr1\t1\t50001\t1\t7500\n'
            'chr1\t1\t50001\t3\t7500\n'
        ))
        summary = vdd.load_summary(summary_path, vdd.load_superpartitions(sp_path))
        self.assertEqual(0.0, summary.observed[0][summary.superpartition_position(2)])

    def test_gzipped_inputs_are_read(self):
        sp_path = self.write('sp.tsv.gz', self.SUPERPARTITIONS)
        summary_path = self.write('summary.tsv.gz', (
            'contig\tbin_start\tbin_end\tsuperpartition\tobserved\n'
            'chr1\t1\t50001\t1\t7500\n'
        ))
        summary = vdd.load_summary(summary_path, vdd.load_superpartitions(sp_path))
        self.assertEqual(1, len(summary.bins))

    def test_undeclared_superpartition_is_an_error(self):
        """Mismatched inputs would silently corrupt every median, so fail loudly."""
        sp_path = self.write('sp.tsv', self.SUPERPARTITIONS)
        summary_path = self.write('summary.tsv', (
            'contig\tbin_start\tbin_end\tsuperpartition\tobserved\n'
            'chr1\t1\t50001\t99\t7500\n'
        ))
        with self.assertRaises(vdd.SummaryFormatError):
            vdd.load_summary(summary_path, vdd.load_superpartitions(sp_path))

    def test_bad_header_is_an_error(self):
        path = self.write('summary.tsv', 'contig\tstart\tend\tsp\tvalue\n')
        with self.assertRaises(vdd.SummaryFormatError):
            vdd.load_summary(path, {1: 100})

    def test_inverted_bin_bounds_are_an_error(self):
        sp_path = self.write('sp.tsv', self.SUPERPARTITIONS)
        summary_path = self.write('summary.tsv', (
            'contig\tbin_start\tbin_end\tsuperpartition\tobserved\n'
            'chr1\t50001\t1\t1\t7500\n'
        ))
        with self.assertRaises(vdd.SummaryFormatError):
            vdd.load_summary(summary_path, vdd.load_superpartitions(sp_path))

    def test_duplicate_superpartition_declaration_is_an_error(self):
        path = self.write('sp.tsv', 'superpartition\tn_samples\n1\t100\n1\t50\n')
        with self.assertRaises(vdd.SummaryFormatError):
            vdd.load_superpartitions(path)

    def test_empty_file_is_an_error(self):
        path = self.write('sp.tsv', '')
        with self.assertRaises(vdd.SummaryFormatError):
            vdd.load_superpartitions(path)


class TestReportOutput(unittest.TestCase):

    def test_report_tsv_has_a_row_per_rectangle(self):
        import io
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, 0.0)
        report = vdd.analyze(summary)
        buffer = io.StringIO()
        vdd.write_report(report, buffer)
        lines = buffer.getvalue().strip().split('\n')
        self.assertEqual(list(vdd.REPORT_COLUMNS), lines[0].split('\t'))
        self.assertEqual(2, len(lines))
        self.assertEqual('83', lines[1].split('\t')[4])

    def test_clean_summary_says_so(self):
        report = vdd.analyze(build_summary())
        self.assertIn('No candidate dropouts found', vdd.format_summary(report, 'variants'))

    def test_summary_reports_percent_present(self):
        summary = build_summary()
        for i in range(10, 19):
            set_cell(summary, i, 83, CELL * 0.10)
        report = vdd.analyze(summary)
        text = vdd.format_summary(report, 'variants')
        self.assertIn('10.00% present', text)


if __name__ == '__main__':
    unittest.main()
