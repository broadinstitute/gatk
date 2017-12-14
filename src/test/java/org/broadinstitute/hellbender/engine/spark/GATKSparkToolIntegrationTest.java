package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GATKSparkToolIntegrationTest extends CommandLineProgramTest {

    @CommandLineProgramProperties(
            summary = "Dummy empty command line that requires a reference .",
            oneLineSummary = "empty class",
            programGroup = TestProgramGroup.class
    )
    public static class DummySparkToolRequiresReference extends GATKSparkTool {
        private static final long serialVersionUID = 0L;
        @Override
        public boolean requiresReference() {
            return true;
        }
        @Override
        protected void runTool(JavaSparkContext ctx) {
            //Do-Nothing
        }
    }

    @Override
    public String getTestedClassName() {
        return DummySparkToolRequiresReference.class.getSimpleName();
    }

    @Test(expectedExceptions = UserException.MissingReferenceDictFile.class)
    public void testMissingReferenceDictFileCatch() throws IOException {
        File outFile = GATKBaseTest.createTempFile("bqsrSparkPipelineTest", ".bam");
        final List<String> args = new ArrayList<>();
        args.add("-R");
        args.add("src/test/resources/org/broadinstitute/hellbender/engine/spark/validNoDict.fasta");
        runCommandLine(args);
    }
}
