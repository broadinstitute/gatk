package org.broadinstitute.hellbender.tools.spark.sv.integration;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.FileUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public final class ExtractOriginalAlignmentRecordsByNameSparkIntegrationTest extends CommandLineProgramTest {

    @Test(groups = "sv")
    public void testExtractOriginalAlignmentRecordsByNameSparkRunnableLocal() throws IOException {

        final File tempWorkingDir = BaseTest.createTempDir("extractOriginalAlignmentRecordsByNameSparkIntegrationTest");

        FileUtils.writeLinesToSingleFile(Collections.singleton("asm013903:tig00002").iterator(), tempWorkingDir + "/names.txt");

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
