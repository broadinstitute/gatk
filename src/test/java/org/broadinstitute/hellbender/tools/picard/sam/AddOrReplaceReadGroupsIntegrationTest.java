package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class AddOrReplaceReadGroupsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = getTestDataDir();
    private static final File SAM_FILE = new File(TEST_DATA_DIR, "picard/sam/AddOrReplaceReadGroups/genomic_sorted_5_plus.sam");

    public String getTestedClassName() {
        return AddOrReplaceReadGroups.class.getSimpleName();
    }

    @Test
    public void testAddCommentsToBam() throws Exception {
        final File outputFile = BaseTest.createTempFile("AddOrReplaceReadGroups", ".sam");
        runIt(SAM_FILE, outputFile);
        final File expectedOutBam = new File(TEST_DATA_DIR, "picard/sam/AddOrReplaceReadGroups/genomic_sorted_5_plus.result.sam"); //created using picard 1.130
        SamAssertionUtils.assertSamsEqual(expectedOutBam, outputFile);
    }

    private void runIt(final File inputFile, final File outputFile) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--INPUT", inputFile.getAbsolutePath(),
                "--OUTPUT", outputFile.getAbsolutePath(),
                "--LB", "foo_LB",
                "--PL","foo_PL",
                "--PU", "foo_PU",
                "--SM", "foo_SM",
                "--CN","foo_CN",
                "--DS","foo_DS",
                "--DT","2015-08-21T09:40:00",
                "--PI","157",
                "--PG","foo_PG",
                "--PM","foo_PM"
                )
        );
        runCommandLine(args);
    }

}
