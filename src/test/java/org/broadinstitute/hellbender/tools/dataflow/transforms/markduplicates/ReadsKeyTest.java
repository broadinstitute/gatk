package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;


import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class ReadsKeyTest {
    @DataProvider(name = "sams")
    public Object[][] sams(){
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        final String name = "fred";
        final int refIndex = 0;
        final int alignmentStart = 101;
        final int length = 76;
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, refIndex, alignmentStart, length);
        final boolean reverseStrand = false;
        read.setIsReverseStrand(reverseStrand);
        final String strandStr = reverseStrand ? "r": "f";
        final String key = ReadsKey.FRAGMENT_PREFIX + "-" + "|" + refIndex + "|" + alignmentStart + "|" + strandStr;

        return new Object[][]{
                new Object[]{header, read, key},
        };
    }

    @Test(dataProvider = "sams", groups= {"dataflow"})
    public void countBasesTest(final SAMFileHeader header, final GATKRead read, final String expectedKey){
        final String key = ReadsKey.keyForFragment(header, read);
        Assert.assertEquals(expectedKey, key);
    }}
