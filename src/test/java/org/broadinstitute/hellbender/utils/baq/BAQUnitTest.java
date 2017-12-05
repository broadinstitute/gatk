package org.broadinstitute.hellbender.utils.baq;

import com.google.common.base.Strings;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.Locatable;
import java.nio.file.Path;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public final class BAQUnitTest extends GATKBaseTest {

    private static final int START_CHR = 1;
    private static final int NUM_CHR = 2;
    private static final int CHR_SIZE = 1000;

    public SAMFileHeader createHeader() {
        return ArtificialReadUtils.createArtificialSamHeader(NUM_CHR, START_CHR, CHR_SIZE);
    }

    @Test //regression test for https://github.com/broadinstitute/gatk/issues/1234
    public void testGetReferenceWindowForReadStopsAtContigStart() throws Exception {
        final int bandwidth = 7;
        GATKRead read = ArtificialReadUtils.createArtificialRead("10M");

        final int start = 2;
        final Locatable pos = new SimpleInterval("1", start, start + read.getLength());
        read.setPosition(pos);
        final SimpleInterval referenceWindowForRead = BAQ.getReferenceWindowForRead(read, bandwidth);
        SimpleInterval refWindow = new SimpleInterval("1", 1, read.getEnd() + (bandwidth / 2));  //start is at 1 because we hit the front end of the reference
        Assert.assertEquals(referenceWindowForRead, refWindow);
    }
    private class BAQTest {
        String readBases, refBases;
        byte[] quals, expected;
        String cigar;
        int refOffset;
        int pos;
        ReferenceDataSource rds;

        public BAQTest(String _refBases, String _readBases, String _quals, String _expected) {
            this(_refBases, _readBases, _quals, _expected, null);
        }

        public BAQTest(String _refBases, String _readBases, String _quals, String _expected, ReferenceDataSource rds) {
            this(0, -1, null, _readBases, _refBases, _quals, _expected, rds);
        }

        public BAQTest(int _refOffset, long _pos, String _cigar, String _refBases, String _readBases, String _quals, String _expected, ReferenceDataSource rds) {
            refOffset = _refOffset;
            pos = (int) _pos;
            cigar = _cigar;
            readBases = _readBases;
            refBases = _refBases;
            this.rds = rds;

            quals = new byte[_quals.getBytes().length];
            expected = new byte[_quals.getBytes().length];
            for (int i = 0; i < quals.length; i++) {
                quals[i] = (byte) (_quals.getBytes()[i] - 33);
                expected[i] = (byte) (_expected.getBytes()[i] - 33);
            }
        }

        public String toString() {
            return readBases;
        }

        public GATKRead createRead() {
            GATKRead read = ArtificialReadUtils.createArtificialRead(createHeader(), "foo", 0, pos > 0 ? pos + (refOffset > 0 ? refOffset : 0) : 1, readBases.getBytes(), quals);
            read.setCigar(cigar == null ? String.format("%dM", quals.length) : cigar);
            return read;
        }
    }


    @DataProvider(name = "data")
    public Object[][] createData1() {
        List<BAQTest> params = new ArrayList<>();

        SAMSequenceDictionary dict= new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("1", Integer.MAX_VALUE));

        params.add(new BAQTest(
                "GCTGCTCCTGGTACTGCTGGATGAGGGCCTCGATGAAGCTAAGCTTTTTCTCCTGCTCCTGCGTGATCCGCTGCAG",
                "GCTGCTCCTGGTACTGCTGGATGAGGGCCTCGATGAAGCTAAGCTTTTCCTCCTGCTCCTGCGTGATCCGCTGCAG",
                "?BACCBDDDFFBCFFHHFIHFEIFHIGHHGHBFEIFGIIGEGIIHGGGIHHIIHIIHIIHGICCIGEII@IGIHCG",
                "?BACCBDDDFFBCFFHHFIHFEIFHIGHHGHBFEIFGIIGEGII410..0HIIHIIHIIHGICCIGEII@IGIHCE"));

        params.add(new BAQTest(
                "GCTTTTTCTCCTCCTG",
                "GCTTTTCCTCCTCCTG",
                "IIHGGGIHHIIHHIIH",
                "EI410..0HIIHHIIE"));

        final String refString1 = "AAATTCAAGATTTCAAAGGCTCTTAACTGCTCAAGATAATTTTTTTTTTTTGAGACAGAGTCTTGCTGTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTTGGCTCACTGCAAGCTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCACCCACCACCACGCCTGGCCAATTTTTTTGTATTTTTAGTAGAGATAG";
        final ReferenceDataSource rds1 = new ReferenceMemorySource(new ReferenceBases(refString1.getBytes(), new SimpleInterval("1", 9999807, 10000032)), dict);

        // big and complex, also does a cap from 3 to 4!
        params.add(new BAQTest(-3, 9999810L, "49M1I126M1I20M1I25M",
                refString1,
                "TTCAAGATTTCAAAGGCTCTTAACTGCTCAAGATAATTTTTTTTTTTTGTAGACAGAGTCTTGCTGTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTTGGCTCACTGCAAGCTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCCACCCACCACCACGCCTGGCCTAATTTTTTTGTATTTTTAGTAGAGA",
                ">IHFECEBDBBCBCABABAADBD?AABBACEABABC?>?B>@A@@>A?B3BBC?CBDBAABBBBBAABAABBABDACCCBCDAACBCBABBB:ABDBACBBDCCCCABCDCCBCC@@;?<B@BC;CBBBAB=;A>ACBABBBABBCA@@<?>>AAA<CA@AABBABCC?BB8@<@C<>5;<A5=A;>=64>???B>=6497<<;;<;>2?>BA@??A6<<A59",
                ">EHFECEBDBBCBCABABAADBD?AABBACEABABC?>?B>@A@@>A?838BC?CBDBAABBBBBAABAABBABDACCCBCDAACBCBABBB:ABDBACBBDCCCCABCDCCBCC@@;?<B@BC;CBBBAB=;A>ACBABBBABBCA@@<?>>AAA<CA@AABBABCC?BB8@<@%<>5;<A5=A;>=64>???B;86497<<;;<;>2?>BA@??A6<<A59",
                rds1));

        final String refString2 = "CCGAGTAGCTGGGACTACAGGCACCCACCACCACGCCTGGCC";
        final ReferenceDataSource rds2 = new ReferenceMemorySource(new ReferenceBases(refString2.getBytes(), new SimpleInterval("1", 9999963, 10000004)), dict);

        // now changes
        params.add(new BAQTest(-3, 9999966L, "36M",
                refString2,
                "AGTAGCTGGGACTACAGGCACCCACCACCACGCCTG",
                "A?>>@>AA?@@>A?>A@?>@>>?=>?'>?=>7=?A9",
                "A?>>@>AA?@@>A?>A@?>@>>?=>?'>?=>7=?A9",
                rds2));

        final String refString3 = "CCACCACGCCTGGCCAATTTTTTTGTATTTTTAGTAGAGATA";
        final ReferenceDataSource rds3 = new ReferenceMemorySource(new ReferenceBases(refString3.getBytes(), new SimpleInterval("1", 9999990, 10000031)), dict);

        // raw base qualities are low -- but they shouldn't be capped
        params.add(new BAQTest(-3, 9999993L, "4=13X2=3X1=4X2=4X1=2X",
                refString3,
                "CCACGCTTGGCAAAGTTTTCCGTACGTTTAGCCGAG",
                "33'/(7+270&4),(&&-)$&,%7$',-/61(,6?8",
                "33'/(7+270&4),(&&-)$&,%7$',-/61(,6?8",
                rds3));

        List<Object[]> params2 = new ArrayList<>();
        for (BAQTest x : params) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }


    @Test(dataProvider = "data")
    public void testBAQWithProvidedReference(BAQTest test) {
        if (test.refBases != null) {
            testBAQ(test, false, null);
        }
    }

    @Test(dataProvider = "data")
    public void testBAQWithCigarAndRefLookup(BAQTest test) {
        if (test.cigar != null) {
            testBAQ(test, true, test.rds);
        }
    }

    @Test
    public void testBAQQualRange() {
        BAQ baq = new BAQ(1.0e-3, 0.1, 7, (byte) 4);         // matches current samtools parameters
        final byte ref = (byte) 'A';
        final byte alt = (byte) 'A';

        for (int i = 0; i <= SAMUtils.MAX_PHRED_SCORE; i++) {
            Assert.assertTrue(baq.calcEpsilon(ref, alt, (byte) i) >= 0.0, "Failed to get baq epsilon range");
        }
    }

    @Test
    public void testBAQOverwritesExistingTagWithNull() {
        final Path reference = IOUtils.getPath(hg19_chr1_1M_Reference);
        final ReferenceDataSource rds = new ReferenceFileSource(reference);

        // create a read with a single base off the end of the contig, which cannot be BAQed
        final GATKRead read = ArtificialReadUtils.createArtificialRead(createHeader(), "foo", 0, rds.getSequenceDictionary().getSequence("1").getSequenceLength() + 1, 1);
        read.setBases(new byte[]{(byte) 'A'});
        read.setBaseQualities(new byte[]{(byte) 20});
        read.setCigar("1M");
        read.setAttribute("BQ", "A");

        // try to BAQ and tell it to RECALCULATE AND ADD_TAG
        final BAQ baq = new BAQ(1.0e-3, 0.1, 7, (byte) 4);
        baq.baqRead(read, rds, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);

        // did we remove the existing tag?
        Assert.assertFalse(read.hasAttribute("BQ"));
    }

    private void testBAQ(BAQTest test, boolean lookupWithFasta, ReferenceDataSource rds) {
        final int bandWidth = 7;
        final BAQ baqHMM = new BAQ(1.0e-3, 0.1, bandWidth, (byte) 4);         // matches current samtools parameters

        final GATKRead read = test.createRead();
        final BAQ.BAQCalculationResult result;
        if (lookupWithFasta && test.cigar != null && rds != null) {
            result = baqHMM.calcBAQFromHMM(read, rds);
        } else {
            result = baqHMM.calcBAQFromHMM(read, test.refBases.getBytes(), test.refOffset);
        }
        Assert.assertNotNull(result);

        System.out.println(Strings.repeat("-", 40));
        System.out.println("reads   : " + new String(test.readBases));
        printQuals(System.out, "in-quals:", test.quals, false);
        printQuals(System.out, "bq-quals:", result.bq, false);
        for (int i = 0; i < test.quals.length; i++) {
            Assert.assertTrue(result.bq[i] >= baqHMM.getMinBaseQual() || test.expected[i] < baqHMM.getMinBaseQual(), "BQ < min base quality");
            Assert.assertEquals(result.bq[i], test.expected[i], "Did not see the expected BAQ value at " + i);
        }

    }

    private static void printQuals(PrintStream out, String prefix, byte[] quals, boolean asInt) {
        out.print(prefix);
        for (int i = 0; i < quals.length; i++) {
            if (asInt) {
                out.printf("%2d", (int) quals[i]);
                if (i + 1 != quals.length) out.print(",");
            } else
                out.print((char) (quals[i] + 33));
        }
        out.println();
    }
}
