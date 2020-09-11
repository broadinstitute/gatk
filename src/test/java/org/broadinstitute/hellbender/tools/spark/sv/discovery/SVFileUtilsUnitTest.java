package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class SVFileUtilsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testWriteSAMRecords() {

        final String inputBam = exampleTestDir + "/metrics/exampleMetrics.bam";
        final ReadsDataSource gatkReads = new ReadsDataSource(IOUtils.getPath(inputBam));
        SAMFileHeader expectedHeader = gatkReads.getHeader();
        final List<SAMRecord> expectedSamRecords = new ArrayList<>(100_000);
        gatkReads.forEach(gatkRead -> expectedSamRecords.add(gatkRead.convertToSAMRecord(expectedHeader)));
        final SAMRecordCoordinateComparator samRecordCoordinateComparator = new SAMRecordCoordinateComparator();
        expectedSamRecords.sort(samRecordCoordinateComparator); // sorting again to eliminate the ambiguity about "1-end unmapped pair, which comes first?" issue

        final File tempDirNew = createTempDir("new");
        tempDirNew.deleteOnExit();

        // test output BAM, index file exists, equal contents
        String outputPath = tempDirNew.getAbsolutePath() + "/out.bam";
        SVFileUtils.writeSAMFile(outputPath, expectedSamRecords.iterator(), expectedHeader, true);
        Assert.assertTrue(Files.exists(IOUtils.getPath(outputPath + ".bai"))
                                || Files.exists(IOUtils.getPath(outputPath.replace(".bam", ".bai"))));
        final ReadsDataSource bamOut = new ReadsDataSource(IOUtils.getPath(outputPath));
        final SAMFileHeader bamOutHeader = bamOut.getHeader();
        Assert.assertEquals(bamOutHeader, expectedHeader);
        final List<SAMRecord> bamOutRecords = new ArrayList<>(100_000);
        bamOut.forEach(gatkRead -> bamOutRecords.add(gatkRead.convertToSAMRecord(bamOutHeader)));
        bamOutRecords.sort(samRecordCoordinateComparator);
        Assert.assertEquals(bamOutRecords, expectedSamRecords);

        // test output SAM, equivalent contents
        outputPath = tempDirNew.getAbsolutePath() + "/out.sam";
        SVFileUtils.writeSAMFile(outputPath, expectedSamRecords.iterator(), expectedHeader, false);
        final ReadsDataSource samOut = new ReadsDataSource(IOUtils.getPath(outputPath));
        SAMFileHeader samOutHeader = samOut.getHeader();
        Assert.assertEquals(samOutHeader, expectedHeader);
        final List<SAMRecord> samOutRecords = new ArrayList<>(100_000);
        samOut.forEach(gatkRead -> samOutRecords.add(gatkRead.convertToSAMRecord(samOutHeader)));
        samOutRecords.sort(samRecordCoordinateComparator);
        Assert.assertEquals(expectedSamRecords.size(), samOutRecords.size());
        for (int i = 0; i < expectedSamRecords.size(); ++i) { // direct comparison doesn't apply (BAM vs SAM)
            Assert.assertEquals(samOutRecords.get(i).getSAMString(), expectedSamRecords.get(i).getSAMString());
        }
    }
}
