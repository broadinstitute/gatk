package org.broadinstitute.hellbender.tools.spark.sv.integration;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public final class ExtractOriginalAlignmentRecordsByNameSparkIntegrationTest extends CommandLineProgramTest {

    @Test(groups = "sv")
    public void testExtractOriginalAlignmentRecordsByNameSparkRunnableLocal() throws IOException {

        final File tempWorkingDir = BaseTest.createTempDir("extractOriginalAlignmentRecordsByNameSparkIntegrationTest");

        FileUtils.writeLines(new File(tempWorkingDir, "names.txt"), Collections.singleton("asm013903:tig00002"));

        final SAMFileHeader expectedHeader;
        final List<SAMRecord> expectedRecords;
        try (final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(SVIntegrationTestDataProvider.TEST_CONTIG_SAM))) {
            expectedHeader = readsSource.getHeader();
            expectedRecords = Utils.stream(readsSource.iterator()).filter(r -> r.getName().equals("asm013903:tig00002"))
                    .sorted(Comparator.comparingInt(GATKRead::getAssignedStart)).map(r -> r.convertToSAMRecord(expectedHeader)).collect(Collectors.toList());
        }

        final List<String> args = Arrays.asList("-I", SVIntegrationTestDataProvider.TEST_CONTIG_SAM,
                "-O", tempWorkingDir.getAbsolutePath() + "/names.bam",
                "--readNameFile", tempWorkingDir.getAbsolutePath() + "/names.txt");

        runCommandLine(args);

        try (final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(tempWorkingDir+"/names.bam"))) {
            Assert.assertEquals(readsSource.getHeader(), expectedHeader);
            final List<SAMRecord> samRecords = Utils.stream(readsSource.iterator()).map(r -> r.convertToSAMRecord(expectedHeader)).collect(Collectors.toList());
            Assert.assertEquals(samRecords.stream().map(SAMRecord::getSAMString).collect(Collectors.toList()),
                                expectedRecords.stream().map(SAMRecord::getSAMString).collect(Collectors.toList()));
        }
    }
}
