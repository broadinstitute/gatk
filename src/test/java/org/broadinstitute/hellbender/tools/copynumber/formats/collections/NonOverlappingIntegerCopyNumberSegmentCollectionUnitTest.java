package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Unit test for {@link NonOverlappingIntegerCopyNumberSegmentCollection}.
 */
public final class NonOverlappingIntegerCopyNumberSegmentCollectionUnitTest extends GATKBaseTest {
    private static final int NUM_POINTS_PLACEHOLDER = 2;
    private static final double QUAL_SOME_PLACEHOLDER = 90.0;
    private static final double QUAL_ALL_PLACEHOLDER = 30.0;
    private static final double QUAL_START_PLACEHOLDER = 20.0;
    private static final double QUAL_END_PLACEHOLDER = 25.0;

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalsWithOverlaps() {
        final List<IntegerCopyNumberSegment> intervalsWithOverlaps = Arrays.asList(
                new IntegerCopyNumberSegment(new SimpleInterval("1", 1, 100), new IntegerCopyNumberState(1), new IntegerCopyNumberState(2), NUM_POINTS_PLACEHOLDER, QUAL_SOME_PLACEHOLDER, QUAL_ALL_PLACEHOLDER, QUAL_START_PLACEHOLDER, QUAL_END_PLACEHOLDER),
                new IntegerCopyNumberSegment(new SimpleInterval("1", 100, 200), new IntegerCopyNumberState(1), new IntegerCopyNumberState(2), NUM_POINTS_PLACEHOLDER, QUAL_SOME_PLACEHOLDER, QUAL_ALL_PLACEHOLDER, QUAL_START_PLACEHOLDER, QUAL_END_PLACEHOLDER),
                new IntegerCopyNumberSegment(new SimpleInterval("2", 1, 1), new IntegerCopyNumberState(1), new IntegerCopyNumberState(2), NUM_POINTS_PLACEHOLDER, QUAL_SOME_PLACEHOLDER, QUAL_ALL_PLACEHOLDER, QUAL_START_PLACEHOLDER, QUAL_END_PLACEHOLDER));
        new NonOverlappingIntegerCopyNumberSegmentCollection(AbstractSampleLocatableCollectionUnitTest.METADATA_EXPECTED, intervalsWithOverlaps);
    }
}
