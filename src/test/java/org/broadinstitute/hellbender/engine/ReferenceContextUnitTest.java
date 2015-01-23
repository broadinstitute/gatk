package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class ReferenceContextUnitTest extends BaseTest {

    private static final File TEST_REFERENCE = new File(hg19MiniReference);

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullDataSource() {
        try ( ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE) ) {
            ReferenceContext refContext = new ReferenceContext(null, new GenomeLocParser(reference.getSequenceDictionary()).createGenomeLoc("1", 1, 5));
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullInterval() {
        try ( ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE) ) {
            ReferenceContext refContext = new ReferenceContext(reference, null);
        }
    }

    @DataProvider(name = "WindowlessReferenceIntervalDataProvider")
    public Object[][] getWindowlessReferenceIntervals() {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        GenomeLocParser genomeLocParser = new GenomeLocParser(reference.getSequenceDictionary());
        reference.close();

        return new Object[][] {
                { genomeLocParser.createGenomeLoc("1", 1, 3), "NNN" },
                { genomeLocParser.createGenomeLoc("1", 11041, 11045), "GCAAA" },
                { genomeLocParser.createGenomeLoc("1", 11210, 11220), "CGGTGCTGTGC" },
                { genomeLocParser.createGenomeLoc("2", 9995, 10005), "NNNNNNCGTAT" },
                { genomeLocParser.createGenomeLoc("2", 10001, 10080), "CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCAC" },
                { genomeLocParser.createGenomeLoc("2", 10005, 10084), "TCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCC" },
                { genomeLocParser.createGenomeLoc("2", 15995, 16000), "TGTCAG" }
        };
    }

    @Test(dataProvider = "WindowlessReferenceIntervalDataProvider")
    public void testWindowlessReferenceContext( final GenomeLoc interval, final String expectedBases ) {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        ReferenceContext refContext = new ReferenceContext(reference, interval);

        checkReferenceContextBases(refContext, expectedBases);
        Assert.assertEquals(refContext.getInterval(), interval, "Wrong interval in reference context");
        Assert.assertEquals(refContext.getWindow(), interval, "Window in windowless reference context not equal to original interval");
        Assert.assertEquals(refContext.numWindowLeadingBases(), 0, "Non-zero leading window size in windowless reference context");
        Assert.assertEquals(refContext.numWindowTrailingBases(), 0, "Non-zero trailing window size in windowless reference context");

        reference.close();
    }

    @DataProvider(name = "WindowedReferenceIntervalDataProvider")
    public Object[][] getWindowedReferenceIntervals() {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        GenomeLocParser genomeLocParser = new GenomeLocParser(reference.getSequenceDictionary());
        reference.close();

        return new Object[][] {
                // Window off the start of the contig:
                { genomeLocParser.createGenomeLoc("1", 1, 3), 5, 5, genomeLocParser.createGenomeLoc("1", 1, 8), "NNNNNNNN" },
                // Window in middle of contig with equal, non-zero start and stop offsets
                { genomeLocParser.createGenomeLoc("1", 11041, 11045), 5, 5, genomeLocParser.createGenomeLoc("1", 11036, 11050), "CAGGAGCAAAGTCGC" },
                // Window in middle of contig with start offset only
                { genomeLocParser.createGenomeLoc("1", 11210, 11220), 3, 0, genomeLocParser.createGenomeLoc("1", 11207, 11220), "TCACGGTGCTGTGC" },
                // Window in middle of contig with stop offset only
                { genomeLocParser.createGenomeLoc("2", 9995, 10005), 0, 3, genomeLocParser.createGenomeLoc("2", 9995, 10008), "NNNNNNCGTATCCC" },
                // Window in middle of contig with unequal, non-zero start and stop offsets
                { genomeLocParser.createGenomeLoc("2", 10005, 10084), 3, 8, genomeLocParser.createGenomeLoc("2", 10002, 10092), "GTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCCACACCCAC" },
                // Window off the end of the contig
                { genomeLocParser.createGenomeLoc("2", 15995, 16000), 2, 5, genomeLocParser.createGenomeLoc("2", 15993, 16000), "TGTGTCAG" }
        };
    }

    @Test(dataProvider = "WindowedReferenceIntervalDataProvider")
    public void testWindowedContext( final GenomeLoc interval, final int windowStartOffset, final int windowStopOffset, final GenomeLoc expectedWindow, final String expectedBases ) {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        ReferenceContext refContext = new ReferenceContext(reference, interval, windowStartOffset, windowStopOffset);

        checkReferenceContextBases(refContext, expectedBases);
        Assert.assertEquals(refContext.getInterval(), interval, "Wrong interval in reference context");
        Assert.assertEquals(refContext.getWindow(), expectedWindow, "Window in windowed reference context not equal to expected window");
        Assert.assertEquals(refContext.numWindowLeadingBases(), interval.getStart() - expectedWindow.getStart(),
                            "Leading window size in windowed reference context not equal to expected value");
        Assert.assertEquals(refContext.numWindowTrailingBases(), 0, expectedWindow.getStop() - interval.getStop(),
                            "Trailing window size in windowed reference context not equal to expected value");

        reference.close();
    }

    @Test
    public void testDynamicallyChangingWindow() {
        final ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        final GenomeLocParser parser = new GenomeLocParser(reference.getSequenceDictionary());
        final GenomeLoc interval = parser.createGenomeLoc("1", 11210, 11220);
        final ReferenceContext refContext = new ReferenceContext(reference, interval);
        final String intervalBases = "CGGTGCTGTGC";

        Assert.assertEquals(interval, refContext.getWindow());
        Assert.assertEquals(refContext.numWindowLeadingBases(), 0);
        Assert.assertEquals(refContext.numWindowTrailingBases(), 0);
        checkReferenceContextBases(refContext, intervalBases);

        refContext.setWindow(5, 5);
        Assert.assertEquals(refContext.getWindow(), parser.createGenomeLoc(interval.getContig(), interval.getStart() - 5, interval.getStop() + 5));
        Assert.assertEquals(refContext.numWindowLeadingBases(), 5);
        Assert.assertEquals(refContext.numWindowTrailingBases(), 5);
        checkReferenceContextBases(refContext, "GCTCA" + intervalBases + "CAGGG");

        refContext.setWindow(0, 10);
        Assert.assertEquals(refContext.getWindow(), parser.createGenomeLoc(interval.getContig(), interval.getStart(), interval.getStop() + 10));
        Assert.assertEquals(refContext.numWindowLeadingBases(), 0);
        Assert.assertEquals(refContext.numWindowTrailingBases(), 10);
        checkReferenceContextBases(refContext, intervalBases + "CAGGGCGCCC");

        refContext.setWindow(20, 3);
        Assert.assertEquals(refContext.getWindow(), parser.createGenomeLoc(interval.getContig(), interval.getStart() - 20, interval.getStop() + 3));
        Assert.assertEquals(refContext.numWindowLeadingBases(), 20);
        Assert.assertEquals(refContext.numWindowTrailingBases(), 3);
        checkReferenceContextBases(refContext, "CTACAGGACCCGCTTGCTCA" + intervalBases + "CAG");

        refContext.setWindow(0, 0);
        Assert.assertEquals(interval, refContext.getWindow());
        Assert.assertEquals(refContext.numWindowLeadingBases(), 0);
        Assert.assertEquals(refContext.numWindowTrailingBases(), 0);
        checkReferenceContextBases(refContext, intervalBases);

        reference.close();
    }

    private void checkReferenceContextBases( final ReferenceContext refContext, final String expectedBases ) {
        byte[] contextBases = refContext.getBases();

        List<Byte> contextBasesFromIterator = new ArrayList<Byte>();
        Iterator<Byte> baseIterator = refContext.getBasesIterator();
        while ( baseIterator.hasNext() ) {
            contextBasesFromIterator.add(baseIterator.next());
        }

        Assert.assertEquals(contextBases.length, expectedBases.length(), "Wrong number of bases from refContext.getBases()");
        Assert.assertEquals(contextBasesFromIterator.size(), expectedBases.length(), "Wrong number of bases from refContext.getBasesIterator()");

        byte[] expectedBasesByteArray = expectedBases.getBytes();
        for ( int baseIndex = 0; baseIndex < expectedBases.length(); ++baseIndex ) {
            Assert.assertEquals(contextBases[baseIndex], expectedBasesByteArray[baseIndex], "Base #" + (baseIndex + 1) + " incorrect from refContext.getBases()");
            Assert.assertEquals(contextBasesFromIterator.get(baseIndex).byteValue(), expectedBasesByteArray[baseIndex], "Base #" + (baseIndex + 1) + " incorrect from refContext.getBasesIterator()");
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
        try ( ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE) ) {
            GenomeLoc interval = new GenomeLocParser(reference.getSequenceDictionary()).createGenomeLoc("1", 5, 10);
            ReferenceContext refContext = new ReferenceContext(reference, interval, windowStartOffset, windowStopOffset);
        }
    }

    @Test(dataProvider = "InvalidWindowDataProvider", expectedExceptions = GATKException.class)
    public void testInvalidWindowHandlingPostConstruction( final int windowStartOffset, final int windowStopOffset ) {
        try ( ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE) ) {
            GenomeLoc interval = new GenomeLocParser(reference.getSequenceDictionary()).createGenomeLoc("1", 5, 10);
            ReferenceContext refContext = new ReferenceContext(reference, interval);
            refContext.setWindow(windowStartOffset, windowStopOffset);
        }
    }
}
