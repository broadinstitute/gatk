package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public final class ReferenceContextUnitTest extends GATKBaseTest {

    private static final Path TEST_REFERENCE = IOUtils.getPath(hg19MiniReference);

    @DataProvider(name = "EmptyReferenceContextDataProvider")
    public Object[][] getEmptyReferenceContextData() {
        // Default-constructed ReferenceContexts and ReferenceContexts constructed from null ReferenceDataSources
        // and/or null intervals should behave as empty context objects.
        return new Object[][] {
                { new ReferenceContext() },
                { new ReferenceContext(null, null, 0, 0) },
                { new ReferenceContext(null, new SimpleInterval("1", 1, 1), 0, 0 ) },
                { new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), null) }
        };
    }

    @Test(dataProvider = "EmptyReferenceContextDataProvider")
    public void testEmptyReferenceContext( final ReferenceContext refContext) {
        Assert.assertFalse(refContext.hasBackingDataSource() && refContext.getInterval() != null,
                           "Empty ReferenceContext reports having both a backing data source and an interval");
        Assert.assertEquals(refContext.getBases().length, 0, "Empty ReferenceContext should have returned an empty bases array from getBases()");
        Assert.assertFalse(refContext.iterator().hasNext(), "Empty ReferenceContext should have returned an empty bases iterator from iterator()");
    }

    @DataProvider(name = "WindowlessReferenceIntervalDataProvider")
    public Object[][] getWindowlessReferenceIntervals() {
        return new Object[][] {
                { new SimpleInterval("1", 1, 3), "NNN" },
                { new SimpleInterval("1", 11041, 11045), "GCAAA" },
                { new SimpleInterval("1", 11210, 11220), "CGGTGCTGTGC" },
                { new SimpleInterval("2", 9995, 10005), "NNNNNNCGTAT" },
                { new SimpleInterval("2", 10001, 10080), "CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCAC" },
                { new SimpleInterval("2", 10005, 10084), "TCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCC" },
                { new SimpleInterval("2", 15995, 16000), "TGTCAG" }
        };
    }

    @Test(dataProvider = "WindowlessReferenceIntervalDataProvider")
    public void testWindowlessReferenceContext( final SimpleInterval interval, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            ReferenceContext refContext = new ReferenceContext(reference, interval);

            checkReferenceContextBases(refContext, expectedBases);
            Assert.assertEquals(refContext.getInterval(), interval, "Wrong interval in reference context");
            Assert.assertEquals(refContext.getWindow(), interval, "Window in windowless reference context not equal to original interval");
            Assert.assertEquals(refContext.numWindowLeadingBases(), 0, "Non-zero leading window size in windowless reference context");
            Assert.assertEquals(refContext.numWindowTrailingBases(), 0, "Non-zero trailing window size in windowless reference context");
        }
    }

    @DataProvider(name = "WindowedReferenceIntervalDataProvider")
    public Object[][] getWindowedReferenceIntervals() {
        return new Object[][] {
                // Window off the start of the contig:
                { new SimpleInterval("1", 1, 3), 5, 5, new SimpleInterval("1", 1, 8), "NNNNNNNN" },
                // Window in middle of contig with equal, non-zero start and stop offsets
                { new SimpleInterval("1", 11041, 11045), 5, 5, new SimpleInterval("1", 11036, 11050), "CAGGAGCAAAGTCGC" },
                // Window in middle of contig with start offset only
                { new SimpleInterval("1", 11210, 11220), 3, 0, new SimpleInterval("1", 11207, 11220), "TCACGGTGCTGTGC" },
                // Window in middle of contig with stop offset only
                { new SimpleInterval("2", 9995, 10005), 0, 3, new SimpleInterval("2", 9995, 10008), "NNNNNNCGTATCCC" },
                // Window in middle of contig with unequal, non-zero start and stop offsets
                { new SimpleInterval("2", 10005, 10084), 3, 8, new SimpleInterval("2", 10002, 10092), "GTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCCACACCCAC" },
                // Window off the end of the contig
                { new SimpleInterval("2", 15995, 16000), 2, 5, new SimpleInterval("2", 15993, 16000), "TGTGTCAG" }
        };
    }

    @Test(dataProvider = "WindowedReferenceIntervalDataProvider")
    public void testWindowedContext( final SimpleInterval interval, final int windowStartOffset, final int windowStopOffset, final SimpleInterval expectedWindow, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            ReferenceContext refContext = new ReferenceContext(reference, interval, windowStartOffset, windowStopOffset);

            checkReferenceContextBases(refContext, expectedBases);
            Assert.assertEquals(refContext.getInterval(), interval, "Wrong interval in reference context");
            Assert.assertEquals(refContext.getWindow(), expectedWindow, "Window in windowed reference context not equal to expected window");
            Assert.assertEquals(refContext.numWindowLeadingBases(), interval.getStart() - expectedWindow.getStart(),
                    "Leading window size in windowed reference context not equal to expected value");
            Assert.assertEquals(refContext.numWindowTrailingBases(), 0, expectedWindow.getEnd() - interval.getEnd(),
                    "Trailing window size in windowed reference context not equal to expected value");
        }
    }

    @Test(dataProvider = "WindowedReferenceIntervalDataProvider")
    public void testWindowedContextUsingIntervalObjects( final SimpleInterval interval, final int windowStartOffset, final int windowStopOffset, final SimpleInterval expectedWindow, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            ReferenceContext refContext = new ReferenceContext(reference, interval, expectedWindow);

            checkReferenceContextBases(refContext, expectedBases);
            Assert.assertEquals(refContext.getInterval(), interval, "Wrong interval in reference context");
            Assert.assertEquals(refContext.getWindow(), expectedWindow, "Window in windowed reference context not equal to expected window");
            Assert.assertEquals(refContext.numWindowLeadingBases(), interval.getStart() - expectedWindow.getStart(),
                    "Leading window size in windowed reference context not equal to expected value");
            Assert.assertEquals(refContext.numWindowTrailingBases(), 0, expectedWindow.getEnd() - interval.getEnd(),
                    "Trailing window size in windowed reference context not equal to expected value");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullIntervalAndNonNullWindow() {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            new ReferenceContext(reference, null, new SimpleInterval("1", 1, 3));
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalNotInWindow() {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            new ReferenceContext(reference, new SimpleInterval("1", 1, 3), new SimpleInterval("1", 10, 30));
        }
    }

    @Test
    public void testWindowedContextUsingIntervalObjects_nullWindow() {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            final SimpleInterval ival = new SimpleInterval("1", 1, 3);
            final ReferenceContext refContext = new ReferenceContext(reference, ival, null);
            Assert.assertEquals(refContext.getWindow(), ival);
            Assert.assertEquals(refContext.getInterval(), ival);
        }
    }

    @Test
    public void testDynamicallyChangingWindow() {
        try (final ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            final SimpleInterval interval = new SimpleInterval("1", 11210, 11220);
            final ReferenceContext refContext = new ReferenceContext(reference, interval);
            final String intervalBases = "CGGTGCTGTGC";

            Assert.assertEquals(interval, refContext.getWindow());
            Assert.assertEquals(refContext.numWindowLeadingBases(), 0);
            Assert.assertEquals(refContext.numWindowTrailingBases(), 0);
            checkReferenceContextBases(refContext, intervalBases);
            Assert.assertEquals(refContext.getBase(), intervalBases.getBytes()[0]);
            Assert.assertEquals(refContext.getForwardBases(), intervalBases.getBytes());

            refContext.setWindow(5, 5);
            Assert.assertEquals(refContext.getWindow(), new SimpleInterval(interval.getContig(), interval.getStart() - 5, interval.getEnd() + 5));
            Assert.assertEquals(refContext.numWindowLeadingBases(), 5);
            Assert.assertEquals(refContext.numWindowTrailingBases(), 5);
            checkReferenceContextBases(refContext, "GCTCA" + intervalBases + "CAGGG");
            Assert.assertEquals(refContext.getBase(), intervalBases.getBytes()[0]);
            Assert.assertEquals(refContext.getForwardBases(), (intervalBases+"CAGGG").getBytes());

            refContext.setWindow(0, 10);
            Assert.assertEquals(refContext.getWindow(), new SimpleInterval(interval.getContig(), interval.getStart(), interval.getEnd() + 10));
            Assert.assertEquals(refContext.numWindowLeadingBases(), 0);
            Assert.assertEquals(refContext.numWindowTrailingBases(), 10);
            checkReferenceContextBases(refContext, intervalBases + "CAGGGCGCCC");
            Assert.assertEquals(refContext.getBase(), intervalBases.getBytes()[0]);
            Assert.assertEquals(refContext.getForwardBases(), (intervalBases+"CAGGGCGCCC").getBytes());

            refContext.setWindow(20, 3);
            Assert.assertEquals(refContext.getWindow(), new SimpleInterval(interval.getContig(), interval.getStart() - 20, interval.getEnd() + 3));
            Assert.assertEquals(refContext.numWindowLeadingBases(), 20);
            Assert.assertEquals(refContext.numWindowTrailingBases(), 3);
            checkReferenceContextBases(refContext, "CTACAGGACCCGCTTGCTCA" + intervalBases + "CAG");
            Assert.assertEquals(refContext.getBase(), intervalBases.getBytes()[0]);
            Assert.assertEquals(refContext.getForwardBases(), (intervalBases+"CAG").getBytes());

            refContext.setWindow(0, 0);
            Assert.assertEquals(interval, refContext.getWindow());
            Assert.assertEquals(refContext.numWindowLeadingBases(), 0);
            Assert.assertEquals(refContext.numWindowTrailingBases(), 0);
            checkReferenceContextBases(refContext, intervalBases);
            Assert.assertEquals(refContext.getBase(), intervalBases.getBytes()[0]);
            Assert.assertEquals(refContext.getForwardBases(), intervalBases.getBytes());
        }
    }

    @Test
    public void testGetBasesStaticWindow() {
        try (final ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            final SimpleInterval interval = new SimpleInterval("1", 11210, 11220);
            final ReferenceContext refContext = new ReferenceContext(reference, interval);
            final String intervalBases = "CGGTGCTGTGC";

            checkReferenceContextBasesFromInterval(refContext, intervalBases, interval);
            Assert.assertEquals(refContext.getWindow(), interval);

            checkReferenceContextBasesFromInterval(refContext, "GCTCA" + intervalBases + "CAGGG",
                    new SimpleInterval(interval.getContig(), interval.getStart() - 5, interval.getEnd() + 5)
            );
            Assert.assertEquals(refContext.getWindow(), interval);

            checkReferenceContextBasesFromInterval(refContext, intervalBases + "CAGGGCGCCC",
                    new SimpleInterval(interval.getContig(), interval.getStart() - 0, interval.getEnd() + 10)
            );
            Assert.assertEquals(refContext.getWindow(), interval);

            checkReferenceContextBasesFromInterval(refContext, "CTACAGGACCCGCTTGCTCA" + intervalBases + "CAG",
                    new SimpleInterval(interval.getContig(), interval.getStart() - 20, interval.getEnd() + 3)
            );
            Assert.assertEquals(refContext.getWindow(), interval);
        }
    }

    @DataProvider
    private Object[][] provideForTestCopyConstructor() {
        return new Object[][] {
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("1", 2650, 2650),
                        0,
                        0
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("1", 2650, 2650),
                        3,
                        5
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("1", 2640, 2650),
                        0,
                        0
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("1", 2640, 2650),
                        3,
                        5
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("2", 2650, 2650),
                        3,
                        5
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("2", 2650, 2650),
                        0,
                        0
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("2", 2650, 2660),
                        3,
                        5
                },
                {
                        new SimpleInterval("1", 11210, 11220),
                        new SimpleInterval("2", 2650, 2660),
                        0,
                        0
                },
        };
    }

    @Test(dataProvider = "provideForTestCopyConstructor")
    public void testCopyConstructor(final SimpleInterval originalInterval, final SimpleInterval newInterval, final int leadingBases, final int trailingBases) {
        try (final ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {

            final ReferenceContext refContext = new ReferenceContext(reference, originalInterval, leadingBases, trailingBases);
            Assert.assertEquals(refContext.getInterval(), originalInterval, "Set interval is different from expected interval!");

            final ReferenceContext newRefContext = new ReferenceContext(refContext, newInterval);
            Assert.assertEquals(newRefContext.getInterval(), newInterval, "Set interval is different from expected interval!");

            final SimpleInterval newWindow = newRefContext.getWindow();

            final int newLeadingBases = newInterval.getStart() - newWindow.getStart();
            final int newTrailingBases = newWindow.getEnd() - newInterval.getEnd();

            Assert.assertEquals(newLeadingBases, leadingBases, "New window leading bases are not the same as old window leading bases!");
            Assert.assertEquals(newTrailingBases, trailingBases, "New window trailing bases are not the same as old window trailing bases!");
        }
    }

    private void checkReferenceContextBases( final ReferenceContext refContext, final String expectedBases ) {

        final byte[] contextBases = refContext.getBases();

        final List<Byte> contextBasesFromIterator = new ArrayList<>();
        final Iterator<Byte> baseIterator = refContext.iterator();
        while ( baseIterator.hasNext() ) {
            contextBasesFromIterator.add(baseIterator.next());
        }

        Assert.assertEquals(contextBases.length, expectedBases.length(), "Wrong number of bases from refContext.getBases()");

        final byte[] expectedBasesByteArray = expectedBases.getBytes();
        for ( int baseIndex = 0; baseIndex < expectedBases.length(); ++baseIndex ) {
            Assert.assertEquals(contextBases[baseIndex], expectedBasesByteArray[baseIndex], "Base #" + (baseIndex + 1) + " incorrect from refContext.getBases()");
            Assert.assertEquals(contextBasesFromIterator.get(baseIndex).byteValue(), expectedBasesByteArray[baseIndex], "Base #" + (baseIndex + 1) + " incorrect from refContext.iterator()");
        }
    }

    private void checkReferenceContextBasesFromInterval( final ReferenceContext refContext, final String expectedBases, final SimpleInterval interval ) {

        // Do this once for the interval-based call:
        final byte[] contextBases = refContext.getBases(interval);
        checkReferenceContextBasesFromIntervalHelper(expectedBases, contextBases);

        // Do this again for the leading/trailing bounds-based call:
        final byte[] contextBases2 = refContext.getBases(interval);

        // First check that the two context bases are the same:
        Assert.assertEquals(contextBases2, contextBases);

        // Now check vs the expected values:
        checkReferenceContextBasesFromIntervalHelper(expectedBases, contextBases2);
    }

    private void checkReferenceContextBasesFromIntervalHelper(final String expectedBases, final byte[] contextBases) {
        Assert.assertEquals(contextBases.length, expectedBases.length(), "Wrong number of bases from refContext.getBases()");

        final byte[] expectedBasesByteArray = expectedBases.getBytes();
        for ( int baseIndex = 0; baseIndex < expectedBases.length(); ++baseIndex ) {
            Assert.assertEquals(contextBases[baseIndex], expectedBasesByteArray[baseIndex], "Base #" + (baseIndex + 1) + " incorrect from refContext.getBases()");
        }
    }

    @DataProvider(name = "InvalidWindowDataProvider")
    public Object[][] getInvalidWindows() {
        return new Object[][] {
                // window start offset < 0
                {-1, 1},
                // window stop offset < 0
                {1, -1},
                // window start offset < 0 && window stop offset < 0
                {-1, -1}
        };
    }

    @Test(dataProvider = "InvalidWindowDataProvider", expectedExceptions = GATKException.class)
    public void testInvalidWindowHandlingAtConstruction( final int windowStartOffset, final int windowStopOffset ) {
        try ( ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE) ) {
            SimpleInterval interval = new SimpleInterval("1", 5, 10);
            ReferenceContext refContext = new ReferenceContext(reference, interval, windowStartOffset, windowStopOffset);
        }
    }

    @Test(dataProvider = "InvalidWindowDataProvider", expectedExceptions = GATKException.class)
    public void testInvalidWindowHandlingPostConstruction( final int windowStartOffset, final int windowStopOffset ) {
        try ( ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE) ) {
            SimpleInterval interval = new SimpleInterval("1", 5, 10);
            ReferenceContext refContext = new ReferenceContext(reference, interval);
            refContext.setWindow(windowStartOffset, windowStopOffset);
        }
    }

    @DataProvider(name = "SubintervalDataProvider")
    public Object[][] getSubintervals() {
        return new Object[][] {
                // start (1120x):   01234567890
                // reference bases: CGGTGCTGTGC
                {"1", 11211, 1, "CGG"},
                {"1", 11219, 1, "TGC"},
                {"1", 11217, 2, "CTGTG"}
        };
    }

    @Test(dataProvider = "SubintervalDataProvider")
    public void testGetKmerAround(final String contig, final int start, final int padding, String expectedSubsequence){
        // the interval of a ReferenceContext object is *in*clusive on both ends
        try (final ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
            final SimpleInterval interval = new SimpleInterval(contig, start, start);
            final ReferenceContext refContext = new ReferenceContext(reference, interval);
            final String kmer = refContext.getKmerAround(start, padding);
            Assert.assertEquals(kmer, expectedSubsequence);
        }
    }
}
