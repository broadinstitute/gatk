package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


public final class AssemblyRegionUnitTest extends BaseTest {
    private IndexedFastaSequenceFile seq;
    private String contig;
    private int contigLength;
    private SAMFileHeader header;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg19MiniReference));
        contig = "1";
        contigLength = seq.getSequence(contig).length();
        header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
    }

    @Test
    public void testConstructor(){
        final SimpleInterval loc = new SimpleInterval("1", 10, 20);
        final AssemblyRegion ar = new AssemblyRegion(loc, 2, header);
        Assert.assertEquals(ar.getExtension(), 2);
        Assert.assertEquals(ar.isActive(), true);
        Assert.assertEquals(ar.getSpan(), loc);
        Assert.assertEquals(ar.getHeader(), header);
    }

    @DataProvider(name = "ActionRegionCreationTest")
    public Object[][] makePollingData() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int start : Arrays.asList(1, 10, 100, contigLength - 10, contigLength - 1) ) {
            for ( final int size : Arrays.asList(1, 10, 100, 1000) ) {
                for ( final int ext : Arrays.asList(0, 1, 10, 100) ) {
                    for ( final boolean isActive : Arrays.asList(true, false) ) {
                        for ( final boolean addStates : Arrays.asList(true, false) ) {
                            List<ActivityProfileState> states = null;
                            if ( addStates ) {
                                states = new LinkedList<>();
                                for ( int i = start; i < start + size; i++ ) {
                                    states.add(new ActivityProfileState(new SimpleInterval(contig, i + start, i + start), isActive ? 1.0 : 0.0));
                                }
                            }
                            final SimpleInterval loc = new SimpleInterval(contig, start, start + size - 1);
                            tests.add(new Object[]{loc, states, isActive, ext});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ActionRegionCreationTest")
    public void testCreatingAssemblyRegions(final SimpleInterval loc, final List<ActivityProfileState> supportingStates, final boolean isActive, final int extension) {
        final AssemblyRegion region = new AssemblyRegion(loc, supportingStates, isActive, extension, header);
        Assert.assertFalse(region.isFinalized());
        Assert.assertEquals(region.getSpan(), loc);
        Assert.assertEquals(region.getExtendedSpan().getStart(), Math.max(loc.getStart() - extension, 1));
        Assert.assertEquals(region.getExtendedSpan().getEnd(), Math.min(loc.getEnd() + extension, contigLength));
        Assert.assertEquals(region.getReadSpanLoc().getStart(), Math.max(loc.getStart() - extension, 1));
        Assert.assertEquals(region.getReadSpanLoc().getEnd(), Math.min(loc.getEnd() + extension, contigLength));
        Assert.assertEquals(region.isActive(), isActive);
        Assert.assertEquals(region.getExtension(), extension);
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getSupportingStates(), supportingStates == null ? Collections.emptyList() : supportingStates);
        Assert.assertNotNull(region.toString());

        assertGoodReferenceGetter(region.getAssemblyRegionReference(seq), region.getExtendedSpan(), 0);
        assertGoodReferenceGetter(region.getAssemblyRegionReference(seq, 0), region.getExtendedSpan(), 0);
        assertGoodReferenceGetter(region.getAssemblyRegionReference(seq, 10), region.getExtendedSpan(), 10);
        assertGoodReferenceGetter(region.getFullReference(seq), region.getReadSpanLoc(), 0);
        assertGoodReferenceGetter(region.getFullReference(seq, 0), region.getReadSpanLoc(), 0);
        assertGoodReferenceGetter(region.getFullReference(seq, 10), region.getReadSpanLoc(), 10);


        region.setFinalized(false);
        Assert.assertFalse(region.isFinalized());
        region.setFinalized(true);
        Assert.assertTrue(region.isFinalized());
        region.setFinalized(false);
        Assert.assertFalse(region.isFinalized());
    }

    private void assertGoodReferenceGetter(final byte[] actualBytes, final SimpleInterval span, final int padding) {
        final int expectedStart = Math.max(span.getStart() - padding, 1);
        final int expectedStop = Math.min(span.getEnd() + padding, contigLength);
        final byte[] expectedBytes = seq.getSubsequenceAt(span.getContig(), expectedStart, expectedStop).getBases();
        Assert.assertEquals(actualBytes, expectedBytes);
    }

    @DataProvider(name = "AssemblyRegionReads")
    public Object[][] makeAssemblyRegionReads() {
        final List<Object[]> tests = new ArrayList<>();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        for ( final int start : Arrays.asList(1, 10, 100, contigLength - 10, contigLength - 1) ) {
            for ( final int readStartOffset : Arrays.asList(-100, -10, 0, 10, 100) ) {
                for ( final int readSize : Arrays.asList(10, 100, 1000) ) {
                    final SimpleInterval loc = IntervalUtils.trimIntervalToContig(contig, start, start + 10, header.getSequence(contig).getSequenceLength());

                    final int readStart = Math.max(start + readStartOffset, 1);
                    final int readStop = Math.min(readStart + readSize, contigLength);
                    final int readLength = readStop - readStart + 1;
                    if ( readLength > 0 ) {
                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, readStart, readLength);
                        final SimpleInterval readLoc = new SimpleInterval(read);
                        if ( readLoc.overlaps(loc) ) {
                            tests.add(new Object[]{loc, read});
                        }
                    }
                }
            }
        }

        return tests.subList(2,3).toArray(new Object[][]{});       //HACK!
    }

    @Test(dataProvider = "AssemblyRegionReads")
    public void testAssemblyRegionReads(final SimpleInterval loc, final GATKRead read) throws Exception {
        final SimpleInterval expectedSpan = loc.mergeWithContiguous(new SimpleInterval(read));

        final AssemblyRegion region = new AssemblyRegion(loc, null, true, 0, header);
        final AssemblyRegion region2 = new AssemblyRegion(loc, null, true, 0, header);
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getExtendedSpan(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), loc);
        Assert.assertTrue(region.equalsIgnoreReads(region2));

        region.add(read);
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getExtendedSpan(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), expectedSpan);
        Assert.assertTrue(region.equalsIgnoreReads(region2));

        region.clearReads();
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getExtendedSpan(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), loc);
        Assert.assertTrue(region.equalsIgnoreReads(region2));

        region.addAll(Collections.singleton(read));
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getExtendedSpan(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), expectedSpan);
        Assert.assertTrue(region.equalsIgnoreReads(region2));

        region.removeAll(Collections.<GATKRead>emptySet());
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getExtendedSpan(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), expectedSpan);
        Assert.assertTrue(region.equalsIgnoreReads(region2));

        region.removeAll(Collections.singleton(read));
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getExtendedSpan(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), loc);
        Assert.assertTrue(region.equalsIgnoreReads(region2));

        final GATKRead read2 = read.copy();
        read2.setName(read.getName() + ".clone");

        for ( final GATKRead readToKeep : Arrays.asList(read, read2)) {
            region.addAll(Arrays.asList(read, read2));
            final GATKRead readToDiscard = readToKeep == read ? read2 : read;
            region.removeAll(Collections.singleton(readToDiscard));
            Assert.assertEquals(region.getReads(), Arrays.asList(readToKeep));
            Assert.assertEquals(region.size(), 1);
            Assert.assertEquals(region.getExtendedSpan(), loc);
        }
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure bad inputs are properly detected
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "BadReadsTest")
    public Object[][] makeBadReadsTest() {
        List<Object[]> tests = new ArrayList<>();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        tests.add(new Object[]{
                header,
                ArtificialReadUtils.createArtificialRead(header, "read1", 0, 10, 10),
                ArtificialReadUtils.createArtificialRead(header, "read2", 0, 9, 10)});
        tests.add(new Object[]{
                header,
                ArtificialReadUtils.createArtificialRead(header, "read1", 0, 10, 10),
                ArtificialReadUtils.createArtificialRead(header, "read2", 1, 9, 10)});
        tests.add(new Object[]{
                header,
                ArtificialReadUtils.createArtificialRead(header, "read1", 1, 10, 10),
                ArtificialReadUtils.createArtificialRead(header, "read2", 0, 9, 10)});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BadReadsTest", expectedExceptions = IllegalArgumentException.class)
    public void testBadReads(final SAMFileHeader header, final GATKRead read1, final GATKRead read2) {
        final SimpleInterval loc = new SimpleInterval(read1);
        final AssemblyRegion region = new AssemblyRegion(loc, null, true, 0, header);
        region.add(read1);
        region.add(read2);
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure we can properly cut up an assembly region based on engine intervals
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "SplitAssemblyRegion")
    public Object[][] makeSplitAssemblyRegion() {
        final List<Object[]> tests = new ArrayList<>();

        final SimpleInterval whole_span = new SimpleInterval("1", 1, 500);
        final SimpleInterval gl_before = new SimpleInterval("1", 1, 9);
        final SimpleInterval gl_after = new SimpleInterval("1", 250, 500);
        final SimpleInterval gl_diff_contig = new SimpleInterval("2", 40, 50);

        final int regionStart = 10;
        final int regionStop = 100;
        final SimpleInterval region = new SimpleInterval("1", regionStart, regionStop);

        for ( final SimpleInterval noEffect : Arrays.asList(whole_span) )
            tests.add(new Object[]{
                    region,
                    Arrays.asList(noEffect),
                    Arrays.asList(region)});

        for ( final SimpleInterval noOverlap : Arrays.asList(gl_before, gl_after, gl_diff_contig) )
            tests.add(new Object[]{
                    region,
                    Arrays.asList(noOverlap),
                    Arrays.asList()});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 5, 50)),
                Arrays.asList(new SimpleInterval("1", regionStart, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 50, 200)),
                Arrays.asList(new SimpleInterval("1", 50, regionStop))});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 40, 50)),
                Arrays.asList(new SimpleInterval("1", 40, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 20, 30), new SimpleInterval("1", 40, 50)),
                Arrays.asList(new SimpleInterval("1", 20, 30), new SimpleInterval("1", 40, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 1, 30), new SimpleInterval("1", 40, 50)),
                Arrays.asList(new SimpleInterval("1", regionStart, 30), new SimpleInterval("1", 40, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 1, 30), new SimpleInterval("1", 70, 200)),
                Arrays.asList(new SimpleInterval("1", regionStart, 30), new SimpleInterval("1", 70, regionStop))});

        tests.add(new Object[]{region,
                Arrays.asList(new SimpleInterval("1", 1, 30), new SimpleInterval("1", 40, 50), new SimpleInterval("1", 70, 200)),
                Arrays.asList(new SimpleInterval("1", regionStart, 30), new SimpleInterval("1", 40, 50), new SimpleInterval("1", 70, regionStop))});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SplitAssemblyRegion")
    public void testSplitAssemblyRegion(final SimpleInterval regionLoc, final List<SimpleInterval> intervalLocs, final List<SimpleInterval> expectedRegionLocs) {
        for ( final boolean addSubstates : Arrays.asList(true, false) ) {
            final List<ActivityProfileState> states;
            if ( addSubstates ) {
                states = new LinkedList<>();
                for ( int i = 0; i < regionLoc.size(); i++ )
                    states.add(new ActivityProfileState(new SimpleInterval(regionLoc.getContig(), regionLoc.getStart() + i, regionLoc.getStart() + i), 1.0));
            } else {
                states = null;
            }

            final AssemblyRegion region = new AssemblyRegion(regionLoc, states, true, 0, header);
            final List<AssemblyRegion> regions = region.splitAndTrimToIntervals(new HashSet<>(intervalLocs));

            Assert.assertEquals(regions.size(), expectedRegionLocs.size(), "Wrong number of split locations");
            for ( int i = 0; i < expectedRegionLocs.size(); i++ ) {
                final SimpleInterval expected = expectedRegionLocs.get(i);
                final AssemblyRegion actual = regions.get(i);
                Assert.assertEquals(actual.getSpan(), expected, "Bad region after split");
                Assert.assertEquals(actual.isActive(), region.isActive());
                Assert.assertEquals(actual.getExtension(), region.getExtension());
            }
        }
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure we can properly cut up an assembly region based on engine intervals
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "TrimAssemblyRegionData")
    public Object[][] makeTrimAssemblyRegionData() {
        final List<Object[]> tests = new ArrayList<>();

        // fully enclosed within active region
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 15, 16),
                new SimpleInterval("1", 15, 16), 0});

        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 10, 15),
                new SimpleInterval("1", 10, 15), 0});

        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 15, 20),
                new SimpleInterval("1", 15, 20), 0});

        // needs extra padding on the right
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 15, 25),
                new SimpleInterval("1", 15, 20), 5});

        // needs extra padding on the left
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 5, 15),
                new SimpleInterval("1", 10, 15), 5});

        // needs extra padding on both
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 7, 21),
                new SimpleInterval("1", 10, 20), 3});
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 9, 23),
                new SimpleInterval("1", 10, 20), 3});

        // desired span captures everything, so we're returning everything.  Tests that extension is set correctly
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 1, 50),
                new SimpleInterval("1", 10, 20), 10});

        // At the start of the chromosome, potentially a bit weird
        tests.add(new Object[]{
                new SimpleInterval("1", 1, 10), 10,
                new SimpleInterval("1", 1, 50),
                new SimpleInterval("1", 1, 10), 10});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TrimAssemblyRegionData")
    public void testTrimAssemblyRegion(final SimpleInterval regionLoc, final int extension, final SimpleInterval desiredSpan, final SimpleInterval expectedAssemblyRegion, final int expectedExtension) {
        final AssemblyRegion region = new AssemblyRegion(regionLoc, Collections.<ActivityProfileState>emptyList(), true, extension, header);
        final AssemblyRegion trimmed = region.trim(desiredSpan);
        Assert.assertEquals(trimmed.getSpan(), expectedAssemblyRegion, "Incorrect region");
        Assert.assertEquals(trimmed.getExtension(), expectedExtension, "Incorrect region");
    }

    @Test
    public void testCreateFromReadShard() {
        final File testBam = new File(NA12878_20_21_WGS_bam);
        final File reference = new File(b37_reference_20_21);
        final SimpleInterval shardInterval = new SimpleInterval("20", 10000000, 10001000);
        final SimpleInterval paddedShardInterval = new SimpleInterval(shardInterval.getContig(), shardInterval.getStart() - 100, shardInterval.getEnd() + 100);

        // Traversal settings to match the GATK 3/4 HaplotypeCaller settings when the expected output was generated:
        final int minRegionSize = 50;
        final int maxRegionSize = 300;
        final int regionPadding = 100;
        final double activeProbThreshold = 0.002;
        final int maxProbPropagationDistance = 50;

        // This mock evaluator returns exactly the values that the GATK 4 HaplotypeCallerEngine's isActive()
        // method returns for this region. We don't have direct access to the HaplotypeCallerEngine since
        // we're in public, so we need to use a mock.
        final AssemblyRegionEvaluator mockEvaluator = new AssemblyRegionEvaluator() {
            private final List<SimpleInterval> activeSites = Arrays.asList(
                    new SimpleInterval("20", 9999996, 9999996),
                    new SimpleInterval("20", 9999997, 9999997),
                    new SimpleInterval("20", 10000117, 10000117),
                    new SimpleInterval("20", 10000211, 10000211),
                    new SimpleInterval("20", 10000439, 10000439),
                    new SimpleInterval("20", 10000598, 10000598),
                    new SimpleInterval("20", 10000694, 10000694),
                    new SimpleInterval("20", 10000758, 10000758),
                    new SimpleInterval("20", 10001019, 10001019)
            );

            @Override
            public ActivityProfileState isActive( AlignmentContext locusPileup, ReferenceContext referenceContext, FeatureContext featureContext ) {
                final SimpleInterval pileupInterval = new SimpleInterval(locusPileup);
                return activeSites.contains(pileupInterval) ? new ActivityProfileState(pileupInterval, 1.0) : new ActivityProfileState(pileupInterval, 0.0);
            }
        };

        final List<Pair<SimpleInterval, Boolean>> expectedResults = Arrays.asList(
                // Expected assembly region bounds + flag indicating whether the region is active or not.
                //
                // All but the first and last values EXACTLY match GATK 3.4 ActiveRegion traversal over the same
                // bam file (the first and last values are different because we're artificially cropping the input
                // shard for the purposes of the test, but although they're cropped they are also consistent with
                // GATK 3.4's results).
                Pair.of(new SimpleInterval("20", 9999902, 9999953), false),
                Pair.of(new SimpleInterval("20", 9999954, 10000039), true),
                Pair.of(new SimpleInterval("20", 10000040, 10000079), false),
                Pair.of(new SimpleInterval("20", 10000080, 10000154), true),
                Pair.of(new SimpleInterval("20", 10000155, 10000173), false),
                Pair.of(new SimpleInterval("20", 10000174, 10000248), true),
                Pair.of(new SimpleInterval("20", 10000249, 10000401), false),
                Pair.of(new SimpleInterval("20", 10000402, 10000476), true),
                Pair.of(new SimpleInterval("20", 10000477, 10000560), false),
                Pair.of(new SimpleInterval("20", 10000561, 10000635), true),
                Pair.of(new SimpleInterval("20", 10000636, 10000656), false),
                Pair.of(new SimpleInterval("20", 10000657, 10000795), true),
                Pair.of(new SimpleInterval("20", 10000796, 10000981), false),
                Pair.of(new SimpleInterval("20", 10000982, 10001056), true),
                Pair.of(new SimpleInterval("20", 10001057, 10001100), false)
        );


        try ( final ReadsDataSource readsSource = new ReadsDataSource(testBam);
              final ReferenceDataSource refSource = new ReferenceFileSource(reference) ) {

            // Set the shard's read filter to match the GATK 3/4 HaplotypeCaller settings when the expected output was generated:
            final ReadShard shard = new ReadShard(shardInterval, paddedShardInterval, readsSource);
            shard.setReadFilter(new CountingReadFilter("MAPPING_QUALITY", new MappingQualityReadFilter(20))
                    .and(new CountingReadFilter("MAPPING_QUALITY_AVAILABLE", ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE))
                    .and(new CountingReadFilter("MAPPED", ReadFilterLibrary.MAPPED))
                    .and(new CountingReadFilter("PRIMARY_ALIGNMENT", ReadFilterLibrary.PRIMARY_ALIGNMENT))
                    .and(new CountingReadFilter("NOT_DUPLICATE", ReadFilterLibrary.NOT_DUPLICATE))
                    .and(new CountingReadFilter("PASSES_VENDOR_QUALITY_CHECK", ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK))
                    .and(new CountingReadFilter("GOOD_CIGAR", ReadFilterLibrary.GOOD_CIGAR))
                    .and(new CountingReadFilter("WELLFORMED", new WellformedReadFilter(readsSource.getHeader()))));

            final Iterable<AssemblyRegion> assemblyRegions = AssemblyRegion.createFromReadShard(shard, readsSource.getHeader(), new ReferenceContext(refSource, paddedShardInterval), new FeatureContext(null, paddedShardInterval), mockEvaluator, minRegionSize, maxRegionSize, regionPadding, activeProbThreshold, maxProbPropagationDistance);
            int regionCount = 0;
            for ( final AssemblyRegion region : assemblyRegions ) {
                Assert.assertTrue(regionCount < expectedResults.size(), "Too many regions returned from AssemblyRegion.createFromReadShard()");
                Assert.assertEquals(region.getSpan(), expectedResults.get(regionCount).getLeft(), "Wrong interval for region");
                Assert.assertEquals(region.isActive(), expectedResults.get(regionCount).getRight().booleanValue(), "Region incorrectly marked as active/inactive");
                ++regionCount;
            }

            Assert.assertEquals(regionCount, expectedResults.size(), "Too few regions returned from AssemblyRegion.createFromReadShard()");
        }
    }
}