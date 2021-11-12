package org.broadinstitute.hellbender.engine.spark;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIDHeaderLine;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Set;

public class GATKSparkToolUnitTest extends GATKBaseTest {

    @CommandLineProgramProperties(
            summary = "TestGATKSparkToolWithVariants",
            oneLineSummary = "TestGATKSparkToolWithVariants",
            programGroup = TestProgramGroup.class
    )
    public static class TestGATKSparkToolWithVariants extends GATKSparkTool {
        private static final long serialVersionUID = 0L;

        @Override
        protected void runTool(JavaSparkContext ctx) {
            //Do-Nothing
        }
    }
    @Test
    public void testGetDefaultToolVCFHeaderLines() {
        final TestGATKSparkToolWithVariants tool = new TestGATKSparkToolWithVariants();
        final String[] args = {"--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "true"};
        tool.instanceMain(args);

        Set<VCFHeaderLine> stdHeaderLines = tool.getDefaultToolVCFHeaderLines();
        stdHeaderLines.add(VCFHeader.makeHeaderVersionLine(VCFHeader.DEFAULT_VCF_VERSION));
        VCFHeader hdr = new VCFHeader(stdHeaderLines);

        VCFHeaderLine sourceLine = hdr.getOtherHeaderLineUnique("source");
        Assert.assertEquals(sourceLine.getValue(), tool.getClass().getSimpleName());

        // TODO: this is bogus:
        //  1) this isn't the name of the key anymore...
        //  2) there can be more than one GATKCommandLine header line (or at least it should
        //     be possible to have more than one)
        //  3) this header line is not really guaranteed to be a VCFIDHeaderLine, although it is
        VCFIDHeaderLine commandLine = (VCFIDHeaderLine) hdr.getOtherHeaderLineUnique("GATKCommandLine");
        Assert.assertEquals(commandLine.getID(), tool.getClass().getSimpleName());

        String commandLineString = commandLine.toString();
        assertContains(commandLineString,"CommandLine=");
        assertContains(commandLineString,"Version=");
        assertContains(commandLineString,"Date=");
    }
}
