package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class FlowBasedReadUnitTest extends GATKBaseTest {

    @Test
    void testBAMFormatParsing() throws Exception{
        final String    testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/read/flow/reads/";
        final String inputDir = testResourceDir + "/input/";
        final String outputDir = testResourceDir + "/outputs/";

        final Path inputFile = FileSystems.getDefault().getPath(inputDir, "sample.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputFile.toString()));
        final String flowOrder = "GTAC";
        final FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

        final String tempOutputDir = createTempDir("expected_outputs").toString();

        int curRead = 0;
        final Iterator<SAMRecord> sr;
        for ( sr = reader.iterator(), curRead = 0 ; sr.hasNext(); curRead++) {
            final FlowBasedRead fbr = new FlowBasedRead(sr.next(),flowOrder, 12, fbargs);
            fbr.applyAlignment();
            Assert.assertEquals(fbr.totalKeyBases(), fbr.seqLength());

            try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".key.txt")) {
                fbr.writeKey(fos);
            }


            String expectedFile = outputDir + "sample." + curRead + ".key.txt";
            IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".key.txt"), new File(expectedFile));
            try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".matrix.txt")){
                fbr.writeMatrix(fos);
            }
            expectedFile = outputDir + "sample." + curRead + ".matrix.txt";
            IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".matrix.txt"), new File(expectedFile));
        }
    }

    @Test
    void testBAMFormatParsingWithT0() throws Exception{
        final String    testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/read/flow/reads/";
        final String inputDir = testResourceDir + "/input/";
        final String outputDir = testResourceDir + "/outputs/";

        final Path inputFile = FileSystems.getDefault().getPath(inputDir, "sample.t0.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputFile.toString()));
        final String flowOrder = "TGCA";
        final FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
        fbargs.useT0Tag = true;
        final String tempOutputDir = createTempDir("expected_outputs").toString();

        int curRead = 0;
        final Iterator<SAMRecord> sr;
        for ( sr = reader.iterator(), curRead = 0 ; sr.hasNext(); curRead++) {
            final FlowBasedRead fbr = new FlowBasedRead(sr.next(),flowOrder, 12, fbargs);
            fbr.applyAlignment();
            Assert.assertEquals(fbr.totalKeyBases(), fbr.seqLength());

            try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".key.txt")) {
                fbr.writeKey(fos);
            }


            String expectedFile = outputDir + "sample.t0." + curRead + ".key.txt";
            IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".key.txt"), new File(expectedFile));
            try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".matrix.txt")){
                fbr.writeMatrix(fos);
            }
            expectedFile = outputDir + "sample.t0." + curRead + ".matrix.txt";
            IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".matrix.txt"), new File(expectedFile));
        }
    }


    @Test (dataProvider = "makeReads")
    public void testUncertainFlowTrimming(final GATKRead read, int nTrim, String uncertainFlowBase, final byte[] output, final int start, final int end) {
        FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
        fbargs.flowNumUncertainFlows = nTrim;
        fbargs.flowFirstUncertainFlowBase = uncertainFlowBase;
        GATKRead trimmedRead = FlowBasedRead.hardClipUncertainBases(read, "TAGC", fbargs);
        Assert.assertEquals(trimmedRead.getBases(), output);
        Assert.assertEquals(trimmedRead.getStart(), start);
        Assert.assertEquals(trimmedRead.getEnd(), end);
    }

    private GATKRead makeRead(final byte[] bases, final boolean isReverse) {

        byte[] quals = new byte[bases.length];

        final String cigar = String.format("%dM", bases.length);
        GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, cigar);
        read.setPosition("chr1",100);
        read.setIsReverseStrand(isReverse);

        return read;
    }

    @DataProvider(name="makeReads")
    public Object[][] makeMultipleReads(){
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 4, "T", new byte[] {'G','A'}, 104, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 4, "T", new byte[] {'G','A','T'}, 103, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 4, "T", new byte[] {'T','A','G','C'},100,103});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 1, "T", new byte[] {'A','G','C','G','A'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 2, "T", new byte[] {'C','C','G','A','T'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 10, "T", new byte[] {},0,0});

        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 4, "A", new byte[] {'A','G','C','G','A'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 4, "A", new byte[] {'G','A','T'}, 103, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 4, "A", new byte[] {'T','A','G','C','G'},100,104});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 1, "A", new byte[] {'T','A','G','C','G','A'}, 100, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 2, "A", new byte[] {'C','C','G','A','T'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 10, "A", new byte[] {'T','A','G'},100,102});

        return tests.toArray(new Object[][]{});
    }

    @Test
    public void testFlowBasedReadConstructorEnforcesMatrix() {

        // make a non-flow read
        final GATKRead        read = makeRead("AAAAA".getBytes(), false);

        // try to make it into a flow base read object, should fail
        try {
            new FlowBasedRead(read, FlowBasedRead.DEFAULT_FLOW_ORDER, FlowBasedRead.MAX_CLASS, new FlowBasedArgumentCollection());
        } catch (IllegalStateException e) {
            // this is the positive path, should be getting an exception
            return;
        }

        // if here, failed to get an exception
        Assert.fail("read should have generated an exception");
    }

    @Test
    public void testArtificialFlowBasedReadConstruction() {

        // make a non-flow read
        final GATKRead        read = makeRead("AAAAA".getBytes(), false);

        // convert to a flow based read
        ArtificialReadUtils.makeIntoFlowBased(read);

        // try to make it into a flow base read object, should not fail
        new FlowBasedRead(read, FlowBasedRead.DEFAULT_FLOW_ORDER, FlowBasedRead.MAX_CLASS, new FlowBasedArgumentCollection());
    }


}
