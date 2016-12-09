package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class SATagBuilderUnitTests {

    @Test
    public static void testSetReadsAsSupplementalNo() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        GATKRead primarySup = ArtificialReadUtils.createArtificialRead(header, "read1", 1, 10000, new byte[] {}, new byte[] {}, "4M");
        primarySup.setMappingQuality(100);
        primarySup.setAttribute("NM",20);
        primarySup.setIsReverseStrand(true);
        GATKRead secondarySup = ArtificialReadUtils.createArtificialRead(header, "read1", 2, 10001, new byte[] {}, new byte[] {}, "4S");
        secondarySup.setMappingQuality(200);
        List<GATKRead> sups = new ArrayList<>();
        sups.add(secondarySup);

        SATagBuilder.setReadsAsSupplemental(primarySup,sups);
        Assert.assertEquals(primarySup.getAttributeAsString("SA"), "3,10001,+,4S,200,*;");
        Assert.assertEquals(primarySup.isSupplementaryAlignment(), false);
        Assert.assertEquals(secondarySup.getAttributeAsString("SA"), "2,10000,-,4M,100,20;");
        Assert.assertEquals(secondarySup.isSupplementaryAlignment(), true);

        GATKRead tertiarySup = ArtificialReadUtils.createArtificialRead(header, "read1", 3, 10003, new byte[] {}, new byte[] {}, "4D");
        tertiarySup.setMappingQuality(200);
        primarySup.clearAttribute("SA");
        secondarySup.clearAttribute("SA");
        sups.add(tertiarySup);

        SATagBuilder.setReadsAsSupplemental(primarySup,sups);
        Assert.assertEquals(sups.size(), 2);
        Assert.assertEquals(sups.get(0).getAttributeAsString("SA"), "2,10000,-,4M,100,20;4,10003,+,4D,200,*;");
        Assert.assertTrue(sups.get(1).getAttributeAsString("SA").startsWith("2,10000,-,4M,100,20;"));
        Assert.assertTrue(primarySup.getAttributeAsString("SA").contains("4,10003,+,4D,200,*;"));
        Assert.assertTrue(primarySup.getAttributeAsString("SA").contains("3,10001,+,4S,200,*;"));

        GATKRead unmappedSup = ArtificialReadUtils.createArtificialUnmappedRead(header,new byte[] {}, new byte[] {});
        primarySup.clearAttribute("SA");
        sups.clear();
        sups.add(unmappedSup);

        SATagBuilder.setReadsAsSupplemental(primarySup,sups);
        Assert.assertEquals(primarySup.getAttributeAsString("SA"), "*,0,+,*,0,*;");
        Assert.assertEquals(sups.get(0).getAttributeAsString("SA"), "2,10000,-,4M,100,20;");
    }

    @Test
    public static void testRemoveTagFunctionality() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        GATKRead testRead = ArtificialReadUtils.createArtificialRead(header, "read1", 1, 10002, new byte[] {}, new byte[] {}, "4M");
        final String attributeVal = "2,10000,+,4M,100,20;4,10003,+,4D,200,0;2,10000,-,4M,100,20;";
        testRead.setAttribute("SA", attributeVal);
        SATagBuilder builder = new SATagBuilder(testRead);

        // Testing that unrelated reads are not removed from the list
        GATKRead unRelatedRead = ArtificialReadUtils.createArtificialRead(header, "read3", 3, 1000, new byte[] {}, new byte[] {}, "4D");
        builder.removeTag(unRelatedRead).setSATag();
        Assert.assertEquals(testRead.getAttributeAsString("SA"), attributeVal);

        // Testing that ALL instances of a particular tag get removed
        GATKRead relatedRead = ArtificialReadUtils.createArtificialRead(header, "read1", 1, 10000, new byte[] {}, new byte[] {}, "4M");
        builder.removeTag(relatedRead).setSATag();
        Assert.assertEquals(testRead.getAttributeAsString("SA"),"4,10003,+,4D,200,0;");
    }

    @Test
    public void testAddTag() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        GATKRead read = ArtificialReadUtils.createArtificialRead("40M");
        read.setAttribute("SA","2,500,+,3S2=1X2=2S,60,1;");
        SATagBuilder builder = new SATagBuilder(read);

        GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 1, 10000, new byte[] {}, new byte[] {}, "4M");
        read1.setIsSupplementaryAlignment(true);
        GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read2", 2, 10001, new byte[] {}, new byte[] {}, "4S");
        read2.setIsSupplementaryAlignment(true);
        SATagBuilder read2builder = new SATagBuilder(read2);
        GATKRead read3 = ArtificialReadUtils.createArtificialRead(header, "read3", 3, 10003, new byte[] {}, new byte[] {}, "4D");
        read3.setIsSupplementaryAlignment(false);

        builder.addTag(read1).addTag(read3).addTag(read2builder).setSATag();
        Assert.assertEquals(read.getAttributeAsString("SA"),"4,10003,+,4D,0,*;2,500,+,3S2=1X2=2S,60,1;2,10000,+,4M,0,*;3,10001,+,4S,0,*;");
    }

    @Test
    public void testGetArtificalReadsBasedOnSATag() {
        // setup the record
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("1", 1000));
        header.addSequence(new SAMSequenceRecord("2", 1000));
        GATKRead read = ArtificialReadUtils.createArtificialRead(header, "name", "2", 1, "AAAAAAAAAA".getBytes(), "##########".getBytes());
        read.setIsFirstOfPair();
        read.setCigar(TextCigarCodec.decode("10M"));
        read.setMatePosition("2",1);
        read.setIsSupplementaryAlignment(true);//spec says first 'SA' record will be the primary record

        SATagBuilder builder = new SATagBuilder(read);
        // check no alignments if no SA tag */
        Assert.assertEquals(builder.getArtificialReadsBasedOnSATag(header).size(), 0);

        // Tests that parsing existing tags works properly
        read.setAttribute("SA",
                "2,500,+,3S2=1X2=2S,60,1;" +
                        "1,191,-,8M2S,60,0;");
        builder = new SATagBuilder(read);

        List<GATKRead> suppl = builder.getArtificialReadsBasedOnSATag(header);
        Assert.assertNotNull(suppl);
        Assert.assertEquals(suppl.size(), 2);

        for(final GATKRead other: suppl) {
            Assert.assertFalse(other.isUnmapped());
            Assert.assertTrue(other.isPaired());
            Assert.assertFalse(other.mateIsUnmapped());
            Assert.assertEquals(other.getMateContig(),read.getMateContig());

            Assert.assertEquals(other.getName(),read.getName());
        }

        GATKRead other = suppl.get(0);
        Assert.assertFalse(other.isSupplementaryAlignment());//1st of suppl and 'record' is supplementary
        Assert.assertEquals(other.getContig(),"2");
        Assert.assertEquals(other.getStart(),500);
        Assert.assertEquals(other.getMappingQuality(), 60);
        Assert.assertEquals((int)other.getAttributeAsInteger("NM"),1);
        Assert.assertEquals(other.getCigar().toString(),"3S2=1X2=2S");


        other = suppl.get(1);
        Assert.assertTrue(other.isSupplementaryAlignment());
        Assert.assertEquals(other.getContig(),"1");
        Assert.assertEquals(other.getStart(),191);
        Assert.assertTrue(other.isReverseStrand());
        Assert.assertEquals(other.getMappingQuality(), 60);
        Assert.assertEquals((int)other.getAttributeAsInteger("NM"),0);
        Assert.assertEquals(other.getCigar().toString(),"8M2S");

    }

    @Test(dataProvider = "malformedSATags", expectedExceptions = GATKException.class)
    public void testMalformedSATagCrash(GATKRead read) {
        new SATagBuilder(read);
    }

    @DataProvider(name = "malformedSATags")
    public Object[][] generateMalformedSAReads() {
        GATKRead read1 = ArtificialReadUtils.createArtificialRead("75M");
        read1.setAttribute("SA","2,10000,-,4M,100,20,4;4,10003,+,4D,200,0,4;");
        GATKRead read2 = ArtificialReadUtils.createArtificialRead("75M");
        read2.setAttribute("SA","foo");
        GATKRead read3 = ArtificialReadUtils.createArtificialRead("75M");
        read3.setAttribute("SA","2,20000,-,4M,100,4; ");
        GATKRead read4 = ArtificialReadUtils.createArtificialRead("75M");
        read4.setAttribute("SA","2,30000,-,4M,100,4;;4,10003,+,4D,200,,4;");

        return new Object[][] {{read1},{read2},{read3},{read4}};
    }

    @Test
    public void testEmptyNMTagSupport() {
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        GATKRead read = ArtificialReadUtils.createArtificialRead(header, "name", "2", 1, "AAAAAAAAAA".getBytes(), "##########".getBytes());
        read.setAttribute("SA","1,191,-,8M2S,60,*;");
        List<GATKRead> reads = (new SATagBuilder(read)).getArtificialReadsBasedOnSATag(header);
        Assert.assertFalse(reads.get(0).hasAttribute("NM"));
    }

}
