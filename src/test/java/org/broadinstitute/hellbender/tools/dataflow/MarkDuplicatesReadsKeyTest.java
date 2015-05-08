package org.broadinstitute.hellbender.tools.dataflow;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.dataflow.transforms.MarkDuplicatesReadsKey;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class MarkDuplicatesReadsKeyTest {
    @DataProvider(name = "sams")
    public Object[][] sams(){
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();

        final String name = "fred";
        final int refIndex = 0;
        final int alignmentStart = 101;
        final int length = 76;
        final SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, name, refIndex, alignmentStart, length);
        final boolean negativeStrand = false;
        read.setReadNegativeStrandFlag(negativeStrand);
        final String strandStr = negativeStrand ? "r": "f";
        final String key = "-" + "|" + refIndex + "|" + alignmentStart + "|" + strandStr;

        return new Object[][]{
                new Object[]{read, key},
        };
    }

    @Test(dataProvider = "sams", groups= {"dataflow"})
    public void countBasesTest(final SAMRecord sam, final String expectedKey){
        final Read read = ReadConverter.makeRead(sam);
        final String key = MarkDuplicatesReadsKey.keyForFragment(sam.getHeader(), read);
        Assert.assertEquals(expectedKey, key);
    }}
