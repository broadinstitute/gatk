package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class AddOrReplaceReadGroupsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIR = getTestDataDir() + "/picard/sam/AddOrReplaceReadGroups";

    @DataProvider(name="addOrReplaceReadGroupData")
    public Object[][] addOrReplaceReadGroupData() {
        return new Object[][] {
                { "genomic_sorted_5_plus.sam", null, ".sam", "genomic_sorted_5_plus.result.sam" },
                { "multigroup_valid.bam", null, ".bam", "multigroup_valid.result.bam"},             // result generated with picard 1.140
                { "multigroup_valid.cram", "basic.fasta", ".bam", "multigroup_valid.result.cram"},  // result generated with picard 1.140
                { "multigroup_valid.cram", "basic.fasta", ".cram", "multigroup_valid.result.cram"}  // result generated with picard 1.140
        };
    }

    public String getTestedClassName() {
        return AddOrReplaceReadGroups.class.getSimpleName();
    }

    @Test(dataProvider="addOrReplaceReadGroupData")
    public void testAddOrReplaceReadGroups(final String inputFileName, final String referenceFileName, final String outputExtension, final String expectedFileName) throws Exception {
        final File inputFile = new File(TEST_DATA_DIR, inputFileName);
        final File outputFile = BaseTest.createTempFile("AddOrReplaceReadGroups", outputExtension);
        final File referenceFile = (null == referenceFileName) ? null : new File(TEST_DATA_DIR, referenceFileName);
        runIt(inputFile, referenceFile, outputFile);
        final File expectedOutBam = new File(TEST_DATA_DIR, expectedFileName); //created using picard 1.130
        SamAssertionUtils.assertSamsEqual(outputFile, expectedOutBam, ValidationStringency.SILENT, referenceFile);
    }

    private void runIt(final File inputFile, final File referenceFile, final File outputFile) {
        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--input"); args.add(inputFile.getAbsolutePath());
        args.add("--output"); args.add(outputFile.getAbsolutePath());
        if (null != referenceFile) {
            args.add("--R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.add("--LB"); args.add("foo_LB");
        args.add("--PL"); args.add("foo_PL");
        args.add("--PU"); args.add("foo_PU");
        args.add("--SM"); args.add("foo_SM");
        args.add("--CN"); args.add("foo_CN");
        args.add("--DS"); args.add("foo_DS");
        args.add("--DT"); args.add("2015-08-21T09:40:00");
        args.add("--PI"); args.add("157");
        args.add("--PG"); args.add("foo_PG");
        args.add("--PM"); args.add("foo_PM");

        runCommandLine(args.getArgsList());
    }

}
