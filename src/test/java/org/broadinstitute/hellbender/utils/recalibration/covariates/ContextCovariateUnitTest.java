package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

public final class ContextCovariateUnitTest extends GATKBaseTest {
    private ContextCovariate covariate;
    private RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new ContextCovariate(RAC);
    }

    @Test
    public void testSimpleContexts() {
        final Random rnd = Utils.getRandomGenerator();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        for(int i = 0; i < 10; i++) {
            final GATKRead read = ArtificialReadUtils.createRandomRead(header, 1000,false);

            read.setIsReverseStrand(rnd.nextBoolean());
            final GATKRead clippedRead = ReadClipper.clipLowQualEnds(read, RAC.LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
            final ReadCovariates readCovariates = new ReadCovariates(read.getLength(), 1, new CovariateKeyCache());
            covariate.recordValues(read, header, readCovariates, true);

            verifyCovariateArray(readCovariates.getMismatchesKeySet(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
            verifyCovariateArray(readCovariates.getInsertionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
            verifyCovariateArray(readCovariates.getDeletionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
        }
    }


    private static void verifyCovariateArray(final int[][] values, final int contextSize, final GATKRead read, final Covariate contextCovariate, final byte lowQualTail) {
        for (int i = 0; i < values.length; i++) {
            Assert.assertEquals(contextCovariate.formatKey(values[i][0]), expectedContext(read, i, contextSize, lowQualTail), "offset " + i);
        }
    }

    @DataProvider
    Iterator<Object[]> AnnoyingReads() {
        final List<Object[]> tests = new ArrayList<>();
        for (final byte base : "NBACGTUMLA".getBytes()) {
            for (int i = 0; i < 10; i++) {
                tests.add(new Object[]{base, i});
            }
        }
        return tests.iterator();
    }

    @Test(dataProvider = "AnnoyingReads")
    public void testContextsAnnoyingReads(final byte base, final int i) {
        final Random rnd = Utils.getRandomGenerator();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 10,false);
        read.getBasesNoCopy()[i] = base;

        read.setIsReverseStrand(rnd.nextBoolean());
        final GATKRead clippedRead = ReadClipper.clipLowQualEnds(read, RAC.LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
        final ReadCovariates readCovariates = new ReadCovariates(read.getLength(), 1, new CovariateKeyCache());
        covariate.recordValues(read, header, readCovariates, true);

        verifyCovariateArray(readCovariates.getMismatchesKeySet(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
        verifyCovariateArray(readCovariates.getInsertionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
        verifyCovariateArray(readCovariates.getDeletionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
    }


    @DataProvider(name="strandedBytes")
    public Object[][] strandedBytes() {
        return new Object[][]{
                {"AAAAA", new byte[]{10, 11, 12, 11, 10}, "5M", false, 10, "NAAAN"},
                {"AAAAA", new byte[]{10, 11, 12, 11, 10}, "5M", false, 11, "NNANN"},
                {"AAAAA", new byte[]{10, 11, 12, 11, 10}, "5M", false, 12, ""},
                {"TCGAT", new byte[]{10, 11, 12, 11, 10}, "5M", false, 10, "NCGAN"},
                {"TCGAT", new byte[]{10, 11, 12, 11, 10}, "5M", true,  10, "NTCGN"},
        };
    }
    @Test(dataProvider = "strandedBytes")
    public void testStrandedBytes(final String baseStr, final byte[] quals, final String cigar, final boolean neg, final int lowQTail, final String expecteBaseStr){
        final byte[] bases = baseStr.getBytes();

        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, cigar);
        read.setIsReverseStrand(neg);
        final byte[] strandedBaseArray = ContextCovariate.getStrandedClippedBytes(read, (byte)lowQTail);   //note the cast is due to TestNG limitation - can't use byte as type for lowQTail
        final byte[] expected = expecteBaseStr.getBytes();
        Assert.assertEquals(new String(strandedBaseArray), new String(expected));
    }

    @DataProvider(name="strandedOffset")
    public Object[][] strandedOffset() {
        return new Object[][]{
                {false, 10, 20, 10},   //for positive strand offset is the same
                {false, 10, 100, 10},
                {true, 10, 20, 20-10-1},
                {true, 10, 100, 100-10-1},
        };
    }
    @Test(dataProvider = "strandedOffset")
    public void strandedOffset(final boolean isNegativeStrand, final int offset, final int clippedReadLength, final int expectedStrandedOffset){
        final int strandedOffset = ContextCovariate.getStrandedOffset(isNegativeStrand, offset, clippedReadLength);
        Assert.assertEquals(strandedOffset, expectedStrandedOffset);
    }


    static String expectedContext(final GATKRead originalRead, final int offset, final int contextSize, final byte lowQualTail) {
        final byte[] strandedBaseArray = ContextCovariate.getStrandedClippedBytes(originalRead, lowQualTail);
        final int strandedOffset = ContextCovariate.getStrandedOffset(originalRead.isReverseStrand(), offset, strandedBaseArray.length);

        final int offsetOfContextStart = strandedOffset - contextSize + 1;
        if (offsetOfContextStart < 0) {
            return null;
        } else {
            final int offsetOneAfterContextEnd = offsetOfContextStart + contextSize;
            final String strandedBases = new String(strandedBaseArray);
            final String context = strandedBases.substring(offsetOfContextStart, offsetOneAfterContextEnd);

            for (final byte base : context.getBytes()) {
                if (BaseUtils.simpleBaseToBaseIndex(base) == -1) {
                    return null;
                }
            }
            return context;
        }
    }
}
