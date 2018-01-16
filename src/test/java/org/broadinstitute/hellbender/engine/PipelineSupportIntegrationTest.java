package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;

public class PipelineSupportIntegrationTest extends GATKBaseTest{

    private static final String INPUT_BAM = toolsTestDir + "print_reads.sorted.chr1_1.bam";

    @Test
    // Asserting that Picard programs are able to pipe their output without errors
    public void testPipeForPicardTools() {
        File output = createTempFile("testOutput",".bam");
        BaseTest.runProcess(ProcessController.getThreadLocal(), ("./gatk SortSam -I "+INPUT_BAM+" -O /dev/stdout -SO coordinate " +
                "| ./gatk SetNmMdAndUqTags -I /dev/stdin -O "+ output.getAbsolutePath()+ " --CREATE_INDEX true -R "+hg19_chr1_1M_Reference).split(" "));
    }
}