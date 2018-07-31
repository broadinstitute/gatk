package org.broadinstitute.hellbender.engine;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class AssemblyRegionIteratorUnitTest extends GATKBaseTest {

    @DataProvider
    public Object[][] testCorrectRegionsHaveCorrectReadsAndSizeData() {
        return new Object[][] {
                // One large interval in the shard
                { NA12878_20_21_WGS_bam, b37_reference_20_21, Arrays.asList(new SimpleInterval("20", 10000000, 10100000)), 50, 300, 100 },
                // Multiple intervals in the shard, same contig, no overlap
                { NA12878_20_21_WGS_bam, b37_reference_20_21, Arrays.asList(new SimpleInterval("20", 10000000, 10010000), new SimpleInterval("20", 10040000, 10050000), new SimpleInterval("20", 10060000, 10070000)), 50, 300, 100 },
                // Multiple intervals in the shard that overlap when padded
                { NA12878_20_21_WGS_bam, b37_reference_20_21, Arrays.asList(new SimpleInterval("20", 10000000, 10020000), new SimpleInterval("20", 10020050, 10030000)), 50, 300, 100 },
                // Multiple intervals in the shard, on multiple contigs
                { NA12878_20_21_WGS_bam, b37_reference_20_21, Arrays.asList(new SimpleInterval("20", 10000000, 10010000), new SimpleInterval("21", 10013000, 10015000)), 50, 300, 100 },
        };
    }

    /*
     * This test checks that over various intervals, all assembly regions created by the AssemblyRegionIterator
     * have the correct reads, that the reads are stored in the correct order, and that region padding and other
     * settings are respected.
     *
     * We determine this by going by to the original BAM and doing a fresh query for each region, to ensure that
     * the query results match the reads actually in the region.
     */
    @Test(dataProvider = "testCorrectRegionsHaveCorrectReadsAndSizeData")
    public void testRegionsHaveCorrectReadsAndSize( final String reads, final String reference, final List<SimpleInterval> shardIntervals, final int minRegionSize, final int maxRegionSize, final int assemblyRegionPadding ) throws IOException {
        try ( final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(reads));
              final ReferenceDataSource refSource = ReferenceDataSource.of(IOUtils.getPath(reference)) ) {
            final SAMSequenceDictionary readsDictionary = readsSource.getSequenceDictionary();
            final MultiIntervalLocalReadShard readShard = new MultiIntervalLocalReadShard(shardIntervals, assemblyRegionPadding, readsSource);
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            final AssemblyRegionEvaluator evaluator = new HaplotypeCallerEngine(hcArgs, false, false, readsSource.getHeader(), new CachingIndexedFastaSequenceFile(IOUtils.getPath(b37_reference_20_21)), new VariantAnnotatorEngine(new ArrayList<>(), hcArgs.dbsnp.dbsnp, hcArgs.comps, false));
            final ReadCoordinateComparator readComparator = new ReadCoordinateComparator(readsSource.getHeader());

            final List<ReadFilter> readFilters = new ArrayList<>(2);
            readFilters.add(new WellformedReadFilter());
            readFilters.add(new ReadFilterLibrary.MappedReadFilter());
            final CountingReadFilter combinedReadFilter = CountingReadFilter.fromList(readFilters, readsSource.getHeader());
            readShard.setReadFilter(combinedReadFilter);

            final AssemblyRegionIterator iter = new AssemblyRegionIterator(readShard, readsSource.getHeader(), refSource, null, evaluator, minRegionSize, maxRegionSize, assemblyRegionPadding, 0.002, 50, true);

            AssemblyRegion previousRegion = null;
            while ( iter.hasNext() ) {
                final AssemblyRegion region = iter.next();

                Assert.assertTrue(region.getSpan().size() <= maxRegionSize, "region size " + region.getSpan().size() + " exceeds the configured maximum: " + maxRegionSize);

                final int regionContigLength = readsDictionary.getSequence(region.getSpan().getContig()).getSequenceLength();
                final int expectedLeftRegionPadding = region.getSpan().getStart() - assemblyRegionPadding > 0 ? assemblyRegionPadding : region.getSpan().getStart() - 1;
                final int expectedRightRegionPadding = region.getSpan().getEnd() + assemblyRegionPadding <= regionContigLength ? assemblyRegionPadding : regionContigLength - region.getSpan().getEnd();
                Assert.assertEquals(region.getSpan().getStart() - region.getExtendedSpan().getStart(), expectedLeftRegionPadding, "Wrong amount of padding on the left side of the region");
                Assert.assertEquals(region.getExtendedSpan().getEnd() - region.getSpan().getEnd(), expectedRightRegionPadding, "Wrong amount of padding on the right side of the region");
                final SimpleInterval regionInterval = region.getExtendedSpan();
                final List<GATKRead> regionActualReads = region.getReads();

                if ( previousRegion != null ) {
                    Assert.assertTrue(IntervalUtils.isBefore(previousRegion.getSpan(), region.getSpan(), readsDictionary), "Previous assembly region's span is not before the current assembly region's span");
                    Assert.assertEquals(previousRegion.getSpan().getEnd(), region.getSpan().getStart() - 1, "previous and current regions are not contiguous");
                }

                GATKRead previousRead = null;
                for ( final GATKRead currentRead : regionActualReads ) {
                    if ( previousRead != null ) {
                        Assert.assertTrue(readComparator.compare(previousRead, currentRead) <= 0, "Reads are out of order within the assembly region");
                    }
                    previousRead = currentRead;
                }

                try ( final ReadsDataSource innerReadsSource = new ReadsDataSource(IOUtils.getPath(reads)) ) {
                    final List<GATKRead> regionExpectedReads = Lists.newArrayList(innerReadsSource.query(regionInterval)).stream().filter(combinedReadFilter).collect(Collectors.toList());

                    final List<GATKRead> actualNotInExpected = new ArrayList<>();
                    final List<GATKRead> expectedNotInActual = new ArrayList<>();
                    for ( final GATKRead expectedRead : regionExpectedReads ) {
                        if ( ! regionActualReads.contains(expectedRead) ) {
                            expectedNotInActual.add(expectedRead);
                        }
                    }

                    for ( final GATKRead actualRead : regionActualReads ) {
                        if ( ! regionExpectedReads.contains(actualRead) ) {
                            actualNotInExpected.add(actualRead);
                        }
                    }

                    Assert.assertEquals(regionActualReads.size(), regionExpectedReads.size(), "Wrong number of reads in region " + region + " for extended interval " + regionInterval +
                            ". Expected reads not in actual reads: " + expectedNotInActual + ". Actual reads not in expected reads: " + actualNotInExpected);

                    Assert.assertEquals(regionActualReads, regionExpectedReads, "Wrong reads in region " + region + " for extended interval " + regionInterval +
                            ". Expected reads not in actual reads: " + expectedNotInActual + ". Actual reads not in expected reads: " + actualNotInExpected);
                }
            }
        }
    }

    /**
     * An artificial AssemblyRegionEvaluator used to assert that reads with deletions are or are not present in a pileup
     */
    private static class FakeAssertingAssemblyRegionEvaluator implements AssemblyRegionEvaluator {

        private final SimpleInterval locusWithDeletions;
        private final int expectedNumDeletionsAtLocus;

        public FakeAssertingAssemblyRegionEvaluator(final SimpleInterval locusWithDeletions, final int expectedNumDeletionsAtLocus) {
            this.locusWithDeletions = locusWithDeletions;
            this.expectedNumDeletionsAtLocus = expectedNumDeletionsAtLocus;
        }

        @Override
        public ActivityProfileState isActive(AlignmentContext locusPileup, ReferenceContext referenceContext, FeatureContext featureContext) {
            if ( locusPileup.getLocation().equals(locusWithDeletions) ) {
                int deletionCount = 0;
                for ( final PileupElement pileupElement : locusPileup.getBasePileup() ) {
                    if ( pileupElement.isDeletion() ) {
                        ++deletionCount;
                    }
                }

                Assert.assertEquals(deletionCount, expectedNumDeletionsAtLocus, "Wrong number of deletions in pileup at " + locusPileup.getLocation());
            }

            return new ActivityProfileState(new SimpleInterval(locusPileup), 0.0);
        }
    }

    @DataProvider
    public Object[][] testIncludeReadsWithDeletionsInIsActivePileupsData() {
        return new Object[][] {
                { NA12878_20_21_WGS_bam, b37_reference_20_21, new SimpleInterval("20", 10004770, 10004770), true, 29 },
                { NA12878_20_21_WGS_bam, b37_reference_20_21, new SimpleInterval("20", 10004770, 10004770), false, 0 }
        };
    }

    /*
     * A test to prove that the includeReadsWithDeletionsInIsActivePileups argument to the AssemblyRegionIterator constructor
     * actually causes reads with deletions at a locus to be included (or excluded) from the pileup for that locus sent to
     * the isActive() method of the AssemblyRegionEvaluator. Uses a fake AssemblyRegionEvaluator to check this.
     */
    @Test(dataProvider = "testIncludeReadsWithDeletionsInIsActivePileupsData")
    public void testIncludeReadsWithDeletionsInIsActivePileups(final String reads, final String reference, final SimpleInterval deletionInterval, final boolean includeReadsWithDeletionsInIsActivePileups, final int expectedNumDeletions) {
        try ( final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(reads));
              final ReferenceDataSource refSource = ReferenceDataSource.of(IOUtils.getPath(reference)) ) {
            final SAMSequenceDictionary readsDictionary = readsSource.getSequenceDictionary();
            final SimpleInterval shardInterval = deletionInterval.expandWithinContig(50, readsDictionary);
            final MultiIntervalLocalReadShard readShard = new MultiIntervalLocalReadShard(Arrays.asList(shardInterval), 50, readsSource);

            // Set up our fake AssemblyRegionEvaluator to check that the deletionInterval locus contains
            // expectedNumDeletions reads with deletions in its pileup during the call to isActive()
            final AssemblyRegionEvaluator evaluator = new FakeAssertingAssemblyRegionEvaluator(deletionInterval, expectedNumDeletions);

            final List<ReadFilter> readFilters = new ArrayList<>(2);
            readFilters.add(new WellformedReadFilter());
            readFilters.add(new ReadFilterLibrary.MappedReadFilter());
            final CountingReadFilter combinedReadFilter = CountingReadFilter.fromList(readFilters, readsSource.getHeader());
            readShard.setReadFilter(combinedReadFilter);

            final AssemblyRegionIterator iter = new AssemblyRegionIterator(readShard, readsSource.getHeader(), refSource, null, evaluator, 50, 300, 50, 0.002, 50, includeReadsWithDeletionsInIsActivePileups);

            // Pull from the AssemblyRegionIterator to trigger the call into the FakeAssertingAssemblyRegionEvaluator,
            // which does the actual assert on the pileups passed to isActive()
            while ( iter.hasNext() ) {
                final AssemblyRegion region = iter.next();
            }
        }
    }
}
