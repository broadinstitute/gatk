package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class FragmentUtilsUnitTest extends GATKBaseTest {

    private static final byte HIGH_QUALITY = 30;
    private static final byte OVERLAPPING_QUALITY = 20;

    private GATKRead makeOverlappingRead(final String leftFlank, final int leftQual, final String overlapBases,
                                              final byte[] overlapQuals, final String rightFlank, final int rightQual,
                                              final int alignmentStart) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.addReadGroup(new SAMReadGroupRecord("RG1"));
        final String bases = leftFlank + overlapBases + rightFlank;
        final int readLength = bases.length();
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, alignmentStart, readLength);
        final byte[] leftQuals = Utils.dupBytes((byte) leftQual, leftFlank.length());
        final byte[] rightQuals = Utils.dupBytes((byte) rightQual, rightFlank.length());
        final byte[] quals = Utils.concat(leftQuals, overlapQuals, rightQuals);
        read.setCigar(readLength + "M");
        read.setBases(bases.getBytes());
        read.setBaseQualities(quals);
        read.setReadGroup("RG1");
        read.setMappingQuality(60);
        return read;
    }

    @DataProvider(name = "AdjustFragmentsTest")
    public Object[][] createAdjustFragmentsTest() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
            final byte[] overlappingBaseQuals = new byte[overlapSize];
            for ( int i = 0; i < overlapSize; i++ ) {
                overlappingBaseQuals[i] = HIGH_QUALITY;
            }
            final GATKRead read1  = makeOverlappingRead(leftFlank, HIGH_QUALITY, overlappingBases, overlappingBaseQuals, "", HIGH_QUALITY, 1);
            final GATKRead read2  = makeOverlappingRead("", HIGH_QUALITY, overlappingBases, overlappingBaseQuals, rightFlank, HIGH_QUALITY, leftFlank.length() + 1);
            tests.add(new Object[]{read1, read2, overlapSize});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AdjustFragmentsTest")
    public void testAdjustingTwoReads(final GATKRead read1, final GATKRead read2, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(read1, read2);

        for ( int i = 0; i < read1.getLength() - overlapSize; i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], HIGH_QUALITY);
        }
        for ( int i = read1.getLength() - overlapSize; i < read1.getLength(); i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }

        for ( int i = 0; i < overlapSize; i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }
        for ( int i = overlapSize; i < read2.getLength(); i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], HIGH_QUALITY);
        }
    }
}
