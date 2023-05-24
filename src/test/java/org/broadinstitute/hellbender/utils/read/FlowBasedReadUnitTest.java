package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.Iterator;

public class FlowBasedReadUnitTest extends GATKBaseTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=HaplotypeCallerIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

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

            String expectedFile = outputDir + "sample." + curRead + ".key.txt";
            if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
                Assert.assertEquals(fbr.totalKeyBases(), fbr.seqLength());
                try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".key.txt")) {
                    fbr.writeKey(fos);
                }
                IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".key.txt"), new File(expectedFile));
            } else {
                try (FileWriter fos = new FileWriter(expectedFile)) {
                    fbr.writeKey(fos);
                }
            }
            expectedFile = outputDir + "sample." + curRead + ".matrix.txt";

            if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
                try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".matrix.txt")) {
                    fbr.writeMatrix(fos);
                }
                IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".matrix.txt"), new File(expectedFile));
            } else {
                try (FileWriter fos = new FileWriter(expectedFile)) {
                    fbr.writeMatrix(fos);
                }
            }
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
            String expectedFile = outputDir + "sample.t0." + curRead + ".key.txt";

            if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
                Assert.assertEquals(fbr.totalKeyBases(), fbr.seqLength());
                try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".key.txt")) {
                    fbr.writeKey(fos);
                }
                IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".key.txt"), new File(expectedFile));
            } else {
                try (FileWriter fos = new FileWriter( expectedFile )) {
                    fbr.writeKey(fos);
                }
            }

            expectedFile = outputDir + "sample.t0." + curRead + ".matrix.txt";

            if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
                try (FileWriter fos = new FileWriter(tempOutputDir + "/" + curRead + ".matrix.txt")) {
                    fbr.writeMatrix(fos);
                }
                IntegrationTestSpec.assertEqualTextFiles(new File(tempOutputDir + "/" + curRead + ".matrix.txt"), new File(expectedFile));
            } else {
                try (FileWriter fos = new FileWriter(expectedFile)) {
                    fbr.writeMatrix(fos);
                }

            }
        }
    }



    private GATKRead makeRead(final byte[] bases, final boolean isReverse) {

        byte[] quals = new byte[bases.length];

        final String cigar = String.format("%dM", bases.length);
        GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, cigar);
        read.setPosition("chr1",100);
        read.setIsReverseStrand(isReverse);

        return read;
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

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

}
