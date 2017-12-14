package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Random;

import static org.broadinstitute.hellbender.utils.recalibration.covariates.ContextCovariate.getStrandedClippedBytes;
import static org.broadinstitute.hellbender.utils.recalibration.covariates.ContextCovariate.getStrandedOffset;

public final class ContextCovariateUnitTest extends GATKBaseTest {
    ContextCovariate covariate;
    RecalibrationArgumentCollection RAC;

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
            final GATKRead read = ArtificialReadUtils.createRandomRead(header, 1000);
            read.setIsReverseStrand(rnd.nextBoolean());
            final GATKRead clippedRead = ReadClipper.clipLowQualEnds(read, RAC.LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
            final ReadCovariates readCovariates = new ReadCovariates(read.getLength(), 1, new CovariateKeyCache());
            covariate.recordValues(read, header, readCovariates, true);

            verifyCovariateArray(readCovariates.getMismatchesKeySet(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
            verifyCovariateArray(readCovariates.getInsertionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
            verifyCovariateArray(readCovariates.getDeletionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate, RAC.LOW_QUAL_TAIL);
        }
    }

    public static void verifyCovariateArray(int[][] values, int contextSize, GATKRead read, Covariate contextCovariate, final byte lowQualTail) {
        for (int i = 0; i < values.length; i++) {
            Assert.assertEquals(contextCovariate.formatKey(values[i][0]), expectedContext(read, i, contextSize, lowQualTail), "offset " + i);
        }
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
        final byte[] strandedBaseArray = getStrandedClippedBytes(read, (byte)lowQTail);   //note the cast is due to TestNG limitation - can't use byte as type for lowQTail
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
        final int strandedOffset = getStrandedOffset(isNegativeStrand, offset, clippedReadLength);
        Assert.assertEquals(strandedOffset, expectedStrandedOffset);
    }


    public static String expectedContext (final GATKRead originalRead, final int offset, final int contextSize, final byte lowQualTail) {
        final byte[] strandedBaseArray = getStrandedClippedBytes(originalRead, lowQualTail);
        final int strandedOffset = getStrandedOffset(originalRead.isReverseStrand(), offset, strandedBaseArray.length);

        final int offsetOfContextStart = strandedOffset - contextSize + 1;
        if (offsetOfContextStart < 0) {
            return null;
        } else {
            final int offsetOneAfterContextEnd = offsetOfContextStart + contextSize;
            final String strandedBases = new String(strandedBaseArray);
            final String context = strandedBases.substring(offsetOfContextStart, offsetOneAfterContextEnd);
            if (context.contains("N")) {
                return null;
            } else {
                return context;
            }
        }
    }

}
