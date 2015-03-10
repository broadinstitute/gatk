package org.broadinstitute.hellbender.utils.codecs.hapmap;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Unit tests for the HapMap codec
 */
public class HapMapUnitTest extends BaseTest {
    // our sample hapmap file
    private final static File hapMapFile = new File(publicTestDir + "org/broadinstitute/hellbender/utils/codecs/hapmap/" + "genotypes_chr1_ASW_phase3.3_first500.hapmap");
    private final static String knownLine = "rs2185539 C/T chr1 556738 + ncbi_b36 bbs urn:lsid:bbs.hapmap.org:Protocol:Phase3.r3:1 urn:lsid:bbs.hapmap.org:Assay:Phase3.r3_r" +
            "s2185539:1 urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:4 QC+ CC TC TT CT CC CC CC CC CC CC CC CC CC";
    /**
     * test reading the header off of the file.  We take in the file, read off the first line,
     * close the reader, and then ask the HapMap decoder for the header with a new reader.  These should
     * be equal (i.e. they return the same object).
     */
    @Test
    public void testReadHeader() {
        RawHapMapCodec codec = new RawHapMapCodec();
        final LineIterator reader = getLineIterator();
        try {
            String header = reader.next();
            Assert.assertTrue(header.equals(codec.readActualHeader(getLineIterator())));
        } finally {
            codec.close(reader);
        }
    }

    @Test
    public void testKnownRecordConversion() {
        RawHapMapCodec codec = new RawHapMapCodec();
        RawHapMapFeature feature = codec.decode(knownLine);


        // check that the alleles are right
        Assert.assertEquals(feature.getAlleles().length,2);
        Assert.assertTrue("C".equals(feature.getAlleles()[0]));
        Assert.assertTrue("T".equals(feature.getAlleles()[1]));

        // check the name
        Assert.assertTrue("rs2185539".equals(feature.getName()));

        // check the position
        Assert.assertEquals(feature.getStart(),556738);
        Assert.assertEquals(feature.getEnd(),556738);

        // check the contig
        Assert.assertTrue("chr1".equals(feature.getChr()));
                
        // check the assembly, center, protLSID, assayLSID, panelLSID, and qccode
        Assert.assertTrue("ncbi_b36".equals(feature.getAssembly()));
        Assert.assertTrue("bbs".equals(feature.getCenter()));
        Assert.assertTrue("urn:lsid:bbs.hapmap.org:Protocol:Phase3.r3:1".equals(feature.getProtLSID()));
        Assert.assertTrue("urn:lsid:bbs.hapmap.org:Assay:Phase3.r3_rs2185539:1".equals(feature.getAssayLSID()));
        Assert.assertTrue("urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:4".equals(feature.getPanelLSID()));
        Assert.assertTrue("QC+".equals(feature.getQCCode()));

        // check the strand
        Assert.assertEquals(feature.getStrand(),Strand.POSITIVE);

        // check the genotypes
        int x = 0;
        for (; x < feature.getGenotypes().length; x++) {
            switch (x) {
                case 1: Assert.assertTrue("TC".equals(feature.getGenotypes()[x])); break;
                case 2: Assert.assertTrue("TT".equals(feature.getGenotypes()[x])); break;
                case 3: Assert.assertTrue("CT".equals(feature.getGenotypes()[x])); break;
                default: Assert.assertTrue("CC".equals(feature.getGenotypes()[x])); break;
            }
        }
        // assert that we found the correct number of records
        Assert.assertEquals(x,13);
    }

    @Test
    public void testReadCorrectNumberOfRecords() {
        // setup the record for reading our 500 line file (499 records, 1 header line)
        RawHapMapCodec codec = new RawHapMapCodec();
        final LineIterator reader = getLineIterator();

        int count = 0;
        try {
            codec.readHeader(reader);
            while (reader.hasNext()) {
                codec.decode(reader.next());
                ++count;
            }
        } catch (IOException e) {
            Assert.fail("IOException " + e.getMessage());
        } finally {
            codec.close(reader);
        }
        Assert.assertEquals(count,499);
    }

    @Test
    public void testGetSampleNames() {
        // setup the record for reading our 500 line file (499 records, 1 header line)
        RawHapMapCodec codec = new RawHapMapCodec();
        final LineIterator reader = getLineIterator();

        String line;
        try {
            codec.readHeader(reader);
            line = reader.next();
            RawHapMapFeature feature = codec.decode(line);
            Assert.assertEquals(feature.getSampleIDs().length,87);

        } catch (IOException e) {
            Assert.fail("IOException " + e.getMessage());
        } finally {
            codec.close(reader);
        }
    }


    public LineIterator getLineIterator() {
        try {
            return new LineIteratorImpl(LineReaderUtil.fromBufferedStream(new PositionalBufferedStream(new FileInputStream(hapMapFile))));
        } catch (FileNotFoundException e) {
            Assert.fail("Unable to open hapmap file : " + hapMapFile);
        }
        return null; // for intellij, it doesn't know that assert.fail is fatal
    }

}
