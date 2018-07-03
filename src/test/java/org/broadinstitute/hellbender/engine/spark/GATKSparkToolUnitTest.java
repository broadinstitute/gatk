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
        VCFHeader hdr = new VCFHeader(stdHeaderLines);

        VCFHeaderLine sourceLine = hdr.getOtherHeaderLine("source");
        Assert.assertEquals(sourceLine.getValue(), tool.getClass().getSimpleName());

        VCFIDHeaderLine commandLine = (VCFIDHeaderLine) hdr.getOtherHeaderLine("GATKCommandLine");
        Assert.assertEquals(commandLine.getID(), tool.getClass().getSimpleName());

        String commandLineString = commandLine.toString();
        assertContains(commandLineString,"CommandLine=");
        assertContains(commandLineString,"Version=");
        assertContains(commandLineString,"Date=");
    }
}
