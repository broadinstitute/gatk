package org.broadinstitute.hellbender.utils.codecs.sampileup;

import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.utils.BaseUtils.Base.*;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SAMPileupCodecUnitTest extends GATKBaseTest {

    private final static SAMPileupCodec CODEC = new SAMPileupCodec();

    private final static byte TEST_REFERENCE = A.base;

    // make a line reader for one string
    private LineReader makeReader(final String string) {
        return new LineReader() {

            boolean returned = false;

            @Override
            public String readLine() throws IOException {
                if(returned) {
                    return null;
                } else {
                    returned = true;
                    return string;
                }
            }

            @Override
            public void close() {
            }
        };
    }

    @Test
    public void testCanDecode() {
        final String EXTRA_CHAR = "1";
        for(final String ext: SAMPileupCodec.SAM_PILEUP_FILE_EXTENSIONS) {
            testCanDecodeExtension(ext);
            for (final String bcExt: IOUtil.BLOCK_COMPRESSED_EXTENSIONS) {
                testCanDecodeExtension(ext + bcExt);
            }
        }
    }

    private static void testCanDecodeExtension(final String ext) {
        final String EXTRA_CHAR = "1";
        Assert.assertTrue(CODEC.canDecode("filename." + ext));
        Assert.assertTrue(CODEC.canDecode("filename" + EXTRA_CHAR + "." + ext));
        Assert.assertFalse(CODEC.canDecode("filename." + ext + EXTRA_CHAR));
        Assert.assertFalse(CODEC.canDecode("filename" + ext));
    }

    @DataProvider(name = "stringFeature")
    public Object[][] getStringAndFeature() {
        final StringBuilder bases = new StringBuilder();
        final StringBuilder qualities = new StringBuilder();
        final List<SAMPileupElement> elements = new ArrayList<>();
        for(Object[] obj: getOneSAMPileupElementData()) {
            bases.append(obj[0]);
            qualities.append(obj[1]);
            elements.add((SAMPileupElement) obj[2]);
        }
        final SAMPileupFeature expected = new SAMPileupFeature("seq1", 126, TEST_REFERENCE, elements);
        final String pileupString = String.format("seq1\t126\t%s\t%s\t%s\t%s", (char) TEST_REFERENCE, qualities.length(), bases.toString(), qualities.toString());
        return new Object[][]{{pileupString, expected}};
    }

    @Test(dataProvider = "stringFeature")
    private void testDecode(final String pileupString, final SAMPileupFeature expected) throws Exception {
        final SAMPileupFeature feature = CODEC.decode(pileupString);
        Assert.assertEquals(feature.getContig(), expected.getContig());
        Assert.assertEquals(feature.getStart(), expected.getStart());
        Assert.assertEquals(feature.getEnd(), expected.getEnd());
        Assert.assertEquals(feature.getRef(), expected.getRef());
        Assert.assertEquals(feature.getBasesString(), expected.getBasesString());
        Assert.assertEquals(feature.getQualsString(), expected.getQualsString());
    }

    @Test(dataProvider = "stringFeature")
    private void testDecodeLoc(final String pileupString, final SAMPileupFeature expected) throws Exception {
        final Feature feature = CODEC.decodeLoc(new LineIteratorImpl(makeReader(pileupString)));
        Assert.assertEquals(feature.getContig(), expected.getContig());
        Assert.assertEquals(feature.getStart(), expected.getStart());
        Assert.assertEquals(feature.getEnd(), expected.getEnd());
    }

    @DataProvider(name = "parseOneSAMPileupElement")
    public Object[][] getOneSAMPileupElementData() {
        return new Object[][]{
                // non-reference
                {"A", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"a", ";", new SAMPileupElement(A.base, (byte) 26)},
                {"C", "I", new SAMPileupElement(C.base, (byte) 40)},
                {"c", ";", new SAMPileupElement(C.base, (byte) 26)},
                {"T", "I", new SAMPileupElement(T.base, (byte) 40)},
                {"t", ";", new SAMPileupElement(T.base, (byte) 26)},
                {"G", "I", new SAMPileupElement(G.base, (byte) 40)},
                {"g", ";", new SAMPileupElement(G.base, (byte) 26)},
                {"N", "I", new SAMPileupElement(N.base, (byte) 40)},
                {"n", ";", new SAMPileupElement(N.base, (byte) 26)},
                // reference
                {".", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {",", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},
                // starting of the read
                {"^~A", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"^@a", ";", new SAMPileupElement(A.base, (byte) 26)},
                {"^~.", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {"^@,", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},
                // end of the read
                {"A$", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"a$", ";", new SAMPileupElement(A.base, (byte) 26)},
                {".$", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {",$", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},
                // insertion + base
                {"+11CCCCCCCCCCCA", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"+11CCCCCCCCCCCa", ";", new SAMPileupElement(A.base, (byte) 26)},
                {"+11CCCCCCCCCCC.", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {"+11CCCCCCCCCCC,", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},
                // deletion + base
                {"-11CCCCCCCCCCCA", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"-11CCCCCCCCCCCa", ";", new SAMPileupElement(A.base, (byte) 26)},
                {"-11CCCCCCCCCCC.", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {"-11CCCCCCCCCCC,", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},
                // base + insertion
                {"A+11CCCCCCCCCCC", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"a+11CCCCCCCCCCC", ";", new SAMPileupElement(A.base, (byte) 26)},
                {".+11CCCCCCCCCCC", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {",+11CCCCCCCCCCC", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},
                // base + deletion
                {"A-11CCCCCCCCCCC", "I", new SAMPileupElement(A.base, (byte) 40)},
                {"a-11CCCCCCCCCCC", ";", new SAMPileupElement(A.base, (byte) 26)},
                {".-11CCCCCCCCCCC", "I", new SAMPileupElement(TEST_REFERENCE, (byte) 40)},
                {",-11CCCCCCCCCCC", ";", new SAMPileupElement(TEST_REFERENCE, (byte) 26)},

        };
    }

    @Test(dataProvider = "parseOneSAMPileupElement")
    private void testParseOneBaseAndQuality(final String bases, final String qualities, final SAMPileupElement expected) {
        final List<SAMPileupElement> list = CODEC.parseBasesAndQuals(bases, qualities, TEST_REFERENCE);
        Assert.assertEquals(list.size(), 1, "Broken test because of data provider");
        Assert.assertEquals(list.get(0).getBase(), expected.getBase());
        Assert.assertEquals(list.get(0).getBaseQuality(), expected.getBaseQuality());
    }


}