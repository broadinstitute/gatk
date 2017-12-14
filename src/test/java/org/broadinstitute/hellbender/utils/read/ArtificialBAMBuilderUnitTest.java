package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public final class ArtificialBAMBuilderUnitTest extends GATKBaseTest {

    @DataProvider(name = "ArtificialBAMBuilderUnitTestProvider")
    public Object[][] makeArtificialBAMBuilderUnitTestProvider() {
        final List<Object[]> tests = new LinkedList<>();

        final List<Integer> starts = Arrays.asList(
                1 // very start of the chromosome
        );

        for ( final int readLength : Arrays.asList(10, 20) ) {
            for ( final int skips : Arrays.asList(0, 1, 10) ) {
                for ( final int start : starts ) {
                    for ( final int nSamples : Arrays.asList(1, 2) ) {
                        for ( final int nReadsPerLocus : Arrays.asList(1, 10) ) {
                            for ( final int nLoci : Arrays.asList(10, 100, 1000) ) {
                                final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(nReadsPerLocus, nLoci);
                                bamBuilder.setReadLength(readLength);
                                bamBuilder.setSkipNLoci(skips);
                                bamBuilder.setAlignmentStart(start);
                                bamBuilder.createAndSetHeader(nSamples);
                                tests.add(new Object[]{bamBuilder, readLength, skips, start, nSamples, nReadsPerLocus, nLoci});
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ArtificialBAMBuilderUnitTestProvider")
    public void testBamProvider(final ArtificialBAMBuilder bamBuilder, int readLength, int skips, int start, int nSamples, int nReadsPerLocus, int nLoci) throws IOException {
        Assert.assertEquals(bamBuilder.getReadLength(), readLength);
        Assert.assertEquals(bamBuilder.getSkipNLoci(), skips);
        Assert.assertEquals(bamBuilder.getAlignmentStart(), start);
        Assert.assertEquals(bamBuilder.getAlignmentEnd(), start + nLoci * (skips + 1) + readLength);
        Assert.assertEquals(bamBuilder.getNSamples(), nSamples);
        Assert.assertEquals(bamBuilder.getnReadsPerLocus(), nReadsPerLocus);
        Assert.assertEquals(bamBuilder.getnLoci(), nLoci);

        Assert.assertEquals(bamBuilder.getSamples().size(), bamBuilder.getNSamples());
        Assert.assertNull(bamBuilder.getReference());

        final List<GATKRead> reads = bamBuilder.makeReads();
        Assert.assertEquals(reads.size(), bamBuilder.expectedNumberOfReads());
        for ( final GATKRead read : reads ) {
            assertGoodRead(read.convertToSAMRecord(bamBuilder.getHeader()), bamBuilder);
        }

        final File bam = bamBuilder.makeTemporaryBAMFile();
        final SamReader reader = SamReaderFactory.makeDefault().open(bam);
        Assert.assertTrue(reader.hasIndex());
        final Iterator<SAMRecord> bamIt = reader.iterator();
        int nReadsFromBam = 0;
        int lastStart = -1;
        while ( bamIt.hasNext() ) {
            final SAMRecord read = bamIt.next();
            assertGoodRead(read, bamBuilder);
            nReadsFromBam++;
            Assert.assertTrue(read.getAlignmentStart() >= lastStart);
            lastStart = read.getAlignmentStart();
        }
        Assert.assertEquals(nReadsFromBam, bamBuilder.expectedNumberOfReads());
    }

    private void assertGoodRead(final SAMRecord read, final ArtificialBAMBuilder bamBuilder) {
        Assert.assertEquals(read.getReadLength(), bamBuilder.getReadLength());
        Assert.assertEquals(read.getReadBases().length, bamBuilder.getReadLength());
        Assert.assertEquals(read.getBaseQualities().length, bamBuilder.getReadLength());
        Assert.assertTrue(read.getStart() >= bamBuilder.getAlignmentStart());
        Assert.assertNotNull(read.getReadGroup());
    }
}


