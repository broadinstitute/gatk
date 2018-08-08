package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class SVFileUtilsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testWriteSAMRecords() {

        final String inputBam = largeFileTestDir + "/sv/testSAMWriter_chr20.bam";
        final ReadsDataSource gatkReads = new ReadsDataSource(IOUtils.getPath(inputBam));
        SAMFileHeader header = gatkReads.getHeader();

        final List<SAMRecord> samRecords = new ArrayList<>(100_000);
        gatkReads.forEach(gatkRead -> samRecords.add(gatkRead.convertToSAMRecord(header)));

        final File tempDirNew = BaseTest.createTempDir("new");
        tempDirNew.deleteOnExit();

        // test output BAM, equal contents, index file exists
        String outputPath = tempDirNew.getAbsolutePath() + "/out.bam";
        SVFileUtils.writeSAMFile(outputPath, samRecords.iterator(), header, true);
        Assert.assertTrue(Files.exists(IOUtils.getPath(outputPath + ".bai"))
                                || Files.exists(IOUtils.getPath(outputPath.replace(".bam", ".bai"))));
        final ReadsDataSource bamOut = new ReadsDataSource(IOUtils.getPath(outputPath));
        SAMFileHeader bamOutHeader = bamOut.getHeader();
        final List<SAMRecord> bamOutRecords = new ArrayList<>(100_000);
        bamOut.forEach(gatkRead -> bamOutRecords.add(gatkRead.convertToSAMRecord(bamOutHeader)));
        Assert.assertEquals(samRecords, bamOutRecords);

        // test output SAM, equivalent contents
        outputPath = tempDirNew.getAbsolutePath() + "/out.sam";
        SVFileUtils.writeSAMFile(outputPath, samRecords.iterator(), header, false);

        final ReadsDataSource samOut = new ReadsDataSource(IOUtils.getPath(outputPath));
        SAMFileHeader samOutHeader = samOut.getHeader();
        final List<SAMRecord> samOutRecords = new ArrayList<>(100_000);
        samOut.forEach(gatkRead -> samOutRecords.add(gatkRead.convertToSAMRecord(samOutHeader)));
        Assert.assertEquals(samRecords.size(), samOutRecords.size());
        for (int i = 0; i < samRecords.size(); ++i) { // direct comparison doesn't apply (BAM vs SAM)
            Assert.assertEquals(samRecords.get(i).getSAMString(), samOutRecords.get(i).getSAMString());
        }
    }
}
