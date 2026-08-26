#!/usr/bin/env python3
"""Unit tests for the zone-retry logic in run_in_hail_cluster.

Compute Engine stockouts are a per-zone condition, but Dataproc's default Auto Zone
placement picks one zone and fails rather than moving on, which turns a transient capacity
shortage into a hard workflow failure. These cover the decision logic for retrying in a
sibling zone; the cluster creation itself is not exercised here.

The module imports google.cloud and ijson at top level, neither of which is needed by the
functions under test, so they are stubbed before import.
"""

import sys
import types
import unittest


def _stub(name, **attributes):
    module = types.ModuleType(name)
    for key, value in attributes.items():
        setattr(module, key, value)
    sys.modules[name] = module
    return module


sys.modules.setdefault('ijson', types.ModuleType('ijson'))
if 'google.cloud.dataproc_v1' not in sys.modules:
    sys.modules.setdefault('google', types.ModuleType('google'))
    _dataproc = types.ModuleType('dataproc_v1')
    _stub('google.cloud', dataproc_v1=_dataproc)
    sys.modules['google.cloud.dataproc_v1'] = _dataproc

import run_in_hail_cluster as runner  # noqa: E402


# The message GCP actually returned for VS-1998, kept verbatim so a change in how these
# are recognised is caught against a real example rather than a paraphrase.
REAL_STOCKOUT = (
    "ERROR: (gcloud.dataproc.clusters.create) Operation "
    "[projects/terra-40d6b12d/regions/us-central1/operations/bfa08809] failed: "
    "Error Code: UNAVAILABLE, errorSource: COMPUTE_ENGINE, Error Message: The zone "
    "'projects/terra-40d6b12d/zones/us-central1-f' does not have enough resources "
    "available to fulfill the request.  'NULL:0/NULL:0/NULL:0 (state:STOCKOUT, "
    "sub-state:STOCKOUT, resource type:compute)'."
)


class TestStockoutDetection(unittest.TestCase):

    def test_recognises_the_real_error(self):
        self.assertTrue(runner.looks_like_stockout(REAL_STOCKOUT))

    def test_recognises_resource_pool_exhausted(self):
        self.assertTrue(runner.looks_like_stockout(
            'Error: ZONE_RESOURCE_POOL_EXHAUSTED_WITH_DETAILS'))

    def test_recognises_the_prose_form_alone(self):
        self.assertTrue(runner.looks_like_stockout(
            "The zone 'us-central1-b' does not have enough resources available"))

    def test_ignores_unrelated_failures(self):
        """Retrying elsewhere would not help these, so they must fail fast."""
        for message in (
            'ERROR: (gcloud.dataproc.clusters.create) PERMISSION_DENIED: insufficient rights',
            "ERROR: Quota 'CPUS' exceeded. Limit: 1000.0 in region us-central1.",
            'ERROR: ALREADY_EXISTS: Cluster projects/x/regions/us-central1/clusters/y',
            'Constraint constraints/compute.vmExternalIpAccess violated',
            '',
        ):
            self.assertFalse(runner.looks_like_stockout(message), msg=message)

    def test_bare_unavailable_is_not_enough(self):
        """UNAVAILABLE also covers transient conditions a different zone would not fix."""
        self.assertFalse(runner.looks_like_stockout(
            'Error Code: UNAVAILABLE, errorSource: COMPUTE_ENGINE'))


class TestZoneResolution(unittest.TestCase):

    def test_absent_means_current_behaviour(self):
        """Empty resolves to no zones, so no --zone flag and Dataproc places the cluster."""
        for value in (None, '', '   '):
            self.assertEqual([], runner.resolve_zones(value, 'us-central1', 'proj'))

    def test_explicit_list_is_used_in_order(self):
        self.assertEqual(
            ['us-central1-a', 'us-central1-c', 'us-central1-f'],
            runner.resolve_zones('us-central1-a,us-central1-c,us-central1-f',
                                 'us-central1', 'proj'),
        )

    def test_whitespace_and_empty_entries_tolerated(self):
        self.assertEqual(
            ['us-central1-a', 'us-central1-b'],
            runner.resolve_zones(' us-central1-a , , us-central1-b ', 'us-central1', 'proj'),
        )

    def test_auto_queries_the_region(self):
        calls = []

        def fake_zones_for_region(region, workspace_project):
            calls.append((region, workspace_project))
            return ['us-central1-a', 'us-central1-b']

        original = runner.zones_for_region
        runner.zones_for_region = fake_zones_for_region
        try:
            self.assertEqual(['us-central1-a', 'us-central1-b'],
                             runner.resolve_zones('auto', 'us-central1', 'proj'))
        finally:
            runner.zones_for_region = original
        self.assertEqual([('us-central1', 'proj')], calls)

    def test_auto_is_case_insensitive(self):
        original = runner.zones_for_region
        runner.zones_for_region = lambda region, project: ['z']
        try:
            for value in ('auto', 'AUTO', ' Auto '):
                self.assertEqual(['z'], runner.resolve_zones(value, 'us-central1', 'proj'))
        finally:
            runner.zones_for_region = original

    def test_auto_failing_falls_back_to_dataproc_placement(self):
        """Being unable to list zones is not a reason to fail the workflow."""
        original = runner.zones_for_region
        runner.zones_for_region = lambda region, project: []
        try:
            self.assertEqual([], runner.resolve_zones('auto', 'us-central1', 'proj'))
        finally:
            runner.zones_for_region = original

    def test_a_region_other_than_us_central1_is_not_hardcoded(self):
        self.assertEqual(
            ['europe-west4-a'],
            runner.resolve_zones('europe-west4-a', 'europe-west4', 'proj'),
        )


class TestRetryPolicy(unittest.TestCase):
    """The decision table the retry loop implements, stated directly."""

    @staticmethod
    def should_retry(output, index, zones):
        return index + 1 < len(zones) and runner.looks_like_stockout(output)

    def test_retries_a_stockout_while_zones_remain(self):
        zones = ['a', 'b', 'c']
        self.assertTrue(self.should_retry(REAL_STOCKOUT, 0, zones))
        self.assertTrue(self.should_retry(REAL_STOCKOUT, 1, zones))

    def test_does_not_retry_past_the_last_zone(self):
        self.assertFalse(self.should_retry(REAL_STOCKOUT, 2, ['a', 'b', 'c']))

    def test_does_not_retry_a_non_capacity_failure(self):
        self.assertFalse(self.should_retry('PERMISSION_DENIED', 0, ['a', 'b', 'c']))

    def test_single_attempt_when_no_zones_configured(self):
        """Preserves today's behaviour for callers that pass no zones."""
        self.assertFalse(self.should_retry(REAL_STOCKOUT, 0, [None]))


if __name__ == '__main__':
    unittest.main()
