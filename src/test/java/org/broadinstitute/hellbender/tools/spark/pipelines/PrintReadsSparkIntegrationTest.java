package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.tools.AbstractPrintReadsIntegrationTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

@Test(groups = "spark")
public final class PrintReadsSparkIntegrationTest extends AbstractPrintReadsIntegrationTest {

    @Test()
    public void testUnSortedRoundTrip() throws Exception {
        // This is a technically incorrectly sam with a header indicating that it is coordinate sorted when it is actually
        // queryname sorted. If the ordering is the same after PrintReadsSpark then it means we aren't automatically sorting the output.
        final File inBam = new File(getTestDataDir(), "print_reads.mismatchedHeader.sam");
        try (ReadsDataSource ds = new ReadsDataSource(inBam.toPath())){
            Assert.assertEquals(ds.getHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        }
        final File outBam = GATKBaseTest.createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

}
