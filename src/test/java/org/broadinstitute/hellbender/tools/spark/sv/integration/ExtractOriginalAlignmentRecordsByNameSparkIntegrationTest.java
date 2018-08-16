package org.broadinstitute.hellbender.tools.spark.sv.integration;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public final class ExtractOriginalAlignmentRecordsByNameSparkIntegrationTest extends CommandLineProgramTest {

    @DataProvider
    private Object[][] createTestData() throws IOException {
        final List<Object[]> data = new ArrayList<>(20);

        final File tempWorkingDir = BaseTest.createTempDir("extractOriginalAlignmentRecordsByNameSparkIntegrationTest");

        FileUtils.writeLines(new File(tempWorkingDir, "names.txt"), Collections.singleton("asm013903:tig00002"));

        final SAMFileHeader expectedHeader;
        final List<SAMRecord> expectedRecordsNormalMatch;
        try (final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(SVIntegrationTestDataProvider.TEST_CONTIG_SAM))) {
            expectedHeader = readsSource.getHeader();
            expectedRecordsNormalMatch = Utils.stream(readsSource.iterator()).filter(r -> r.getName().equals("asm013903:tig00002"))
                    .sorted(Comparator.comparingInt(GATKRead::getAssignedStart)).map(r -> r.convertToSAMRecord(expectedHeader)).collect(Collectors.toList());
        }
        final List<String> normalArgs = Arrays.asList("-I", SVIntegrationTestDataProvider.TEST_CONTIG_SAM,
                "-O", tempWorkingDir.getAbsolutePath() + "/names.bam",
                "--read-name-file", tempWorkingDir.getAbsolutePath() + "/names.txt");
        data.add(new Object[]{IOUtils.getPath(tempWorkingDir+"/names.bam"), normalArgs, expectedHeader, expectedRecordsNormalMatch});

        final List<SAMRecord> expectedRecordsInvertedMatch;
        try (final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(SVIntegrationTestDataProvider.TEST_CONTIG_SAM))) {
            expectedRecordsInvertedMatch = Utils.stream(readsSource.iterator()).filter(r -> !r.getName().equals("asm013903:tig00002"))
                    .sorted(Comparator.comparingInt(GATKRead::getAssignedStart)).map(r -> r.convertToSAMRecord(expectedHeader)).collect(Collectors.toList());
        }
        final List<String> invertArgs = Arrays.asList("-I", SVIntegrationTestDataProvider.TEST_CONTIG_SAM,
                "-O", tempWorkingDir.getAbsolutePath() + "/namesInverted.bam",
                "--read-name-file", tempWorkingDir.getAbsolutePath() + "/names.txt",
                "--invert-match");

        data.add(new Object[]{IOUtils.getPath(tempWorkingDir+"/namesInverted.bam"), invertArgs, expectedHeader, expectedRecordsInvertedMatch});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "createTestData")
    public void testExtractOriginalAlignmentRecordsByNameSpark(final Path bamPath, final List<String> args,
                                                               final SAMFileHeader expectedHeader, final List<SAMRecord> expectedRecords) {
        runCommandLine(args);

        try (final ReadsDataSource readsSource = new ReadsDataSource(bamPath)) {
            Assert.assertEquals(readsSource.getHeader(), expectedHeader);
            final List<SAMRecord> samRecords = Utils.stream(readsSource.iterator()).map(r -> r.convertToSAMRecord(expectedHeader)).collect(Collectors.toList());
            Assert.assertEquals(samRecords.stream().map(SAMRecord::getSAMString).collect(Collectors.toList()),
                    expectedRecords.stream().map(SAMRecord::getSAMString).collect(Collectors.toList()));
        }
    }
}
