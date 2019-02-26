package org.broadinstitute.hellbender.tools.walkers.bqsr;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public final class ApplyBQSRIntegrationTest extends AbstractApplyBQSRIntegrationTest {

    @Override
    public String getTestedClassName() {
        return ApplyBQSR.class.getSimpleName();
    }

    @DataProvider(name = "MiniApplyBQSRTest")
    public Object[][] createMiniABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MiniApplyBQSRTest")
    public void testApplyBQSRPath(ABQSRTest params) throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path outPath = jimfs.getPath("applyBQSRTest"+params.outputExtension);

            final ArrayList<String> args = new ArrayList<>();
            Path refPath = null;

            args.add("-I");
            args.add(new File(params.bam).getAbsolutePath());
            args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
            args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
            args.add("-O"); args.add(outPath.toUri().toString());
            if (params.reference != null) {
                File refFile = new File(params.reference);
                args.add("-R"); args.add(refFile.getAbsolutePath());
                refPath = refFile.toPath();
            }
            if (params.args != null) {
                Stream.of(params.args).forEach(arg -> args.add(arg));
            }

            runCommandLine(args);

            SamAssertionUtils.assertSamsEqual(outPath, new File(params.expectedFile).toPath(), refPath);
        }
    }

    @Test(dataProvider = "ApplyBQSRTest", groups={"bucket"})
    public void testApplyBQSRCloud(ABQSRTest params) throws IOException {
        // getTempFilePath also deletes the file on exit.
        final String outString = BucketUtils.getTempFilePath(getGCPTestStaging() + "tmp/testApplyBQSRCloud",  params.outputExtension);
        final Path outPath = BucketUtils.getPathOnGcs(outString);
        final ArrayList<String> args = new ArrayList<>();
        Path refPath = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O");
        args.add(outString);
        if (params.reference != null) {
            File refFile = new File(params.reference);
            args.add("-R");
            args.add(refFile.getAbsolutePath());
            refPath = refFile.toPath();
        }
        if (params.args != null) {
            Stream.of(params.args).forEach(arg -> args.add(arg));
        }

        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outPath, new File(params.expectedFile).toPath(), refPath);
    }

    @Test
    public void testAddingPG() throws IOException {
        final File inFile = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = GATKBaseTest.createTempFile("testAddingPG", ".bam");
        final String[] args = new String[] {
                "--input", inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--use-original-qualities",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //The expected output is actually the same as inputs for this read (this ignores the PGs header)
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has no GATK ApplyBQSR in headers
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));

        //output has a GATK ApplyBQSR in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));
    }
}
