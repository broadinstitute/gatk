package org.broadinstitute.hellbender;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.CRAMIssue8768Detector;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;

public class CRAMIssue8768DetectorIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "cramAnalysisTestCases")
    public Object[][] getFileDiagnosticsTestCases() {
        return new Object[][]{
                // local file tests
                {
                        // test file created by rewriting the htsjdk test file
                        // src/test/resources/htsjdk/samtools/cram/mitoAlignmentStartTest.cram, using a version
                        // of GATK that has bug https://github.com/broadinstitute/gatk/issues/8768. the rewrite causes
                        // the file to have corrupt read bases, so we can use it to test the detector
                        // 1 bad container
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.bug8768.cram",
                        "src/test/resources/filediagnostics/mitoCorrupt.bug8768.fa",
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "false")),
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.false.bug8768.txt",
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.false.bug8768.tsv",
                },
                {
                        // test file created by rewriting the htsjdk test file
                        // src/test/resources/htsjdk/samtools/cram/mitoAlignmentStartTest.cram, using a version
                        // of GATK that has bug https://github.com/broadinstitute/gatk/issues/8768. the rewrite causes
                        // the file to have corrupt read bases, so we can use it to test the detector
                        // 1 bad container
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.bug8768.cram",
                        "src/test/resources/filediagnostics/mitoCorrupt.bug8768.fa",
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "true")),
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.true.bug8768.txt",
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.true.bug8768.tsv",
                },
                {
                        // 0 bad containers
                        NA12878_20_21_WGS_cram,
                        b37_reference_20_21,
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "false")),
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.false.bug8768.txt",
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.false.bug8768.tsv",
                },
                {
                        // 0 bad containers
                        NA12878_20_21_WGS_cram,
                        b37_reference_20_21,
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "true")),
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.true.bug8768.txt",
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.true.bug8768.tsv",
                },
                {
                        // test file created by rewriting the htsjdk test file
                        // src/test/resources/htsjdk/samtools/cram/mitoAlignmentStartTest.cram using a version
                        // of GATK that has bug https://github.com/broadinstitute/gatk/issues/8768, along with code
                        // to force the file to have 2 bad containers by forcing numerous reads to be aligned to
                        // position 1.
                        // 2 bad containers
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.2BadContainers.bug8768.cram",
                        "src/test/resources/filediagnostics/mitoCorrupt.bug8768.fa",
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "false")),
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.2BadContainers.false.bug8768.txt",
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.2BadContainers.false.bug8768.tsv",
                },
                {
                        // test file created by rewriting the htsjdk test file
                        // src/test/resources/htsjdk/samtools/cram/mitoAlignmentStartTest.cram using a version
                        // of GATK that has bug https://github.com/broadinstitute/gatk/issues/8768, along with code
                        // to force the file to have 2 bad containers by forcing numerous reads to be aligned to
                        // position 1.
                        // 3 bad containers
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.3BadContainers.bug8768.cram",
                        "src/test/resources/filediagnostics/mitoCorrupt.bug8768.fa",
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "false")),
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.3BadContainers.false.bug8768.txt",
                        "src/test/resources/filediagnostics/thisFileIsCorruptMito.3BadContainers.false.bug8768.tsv",
                },
                // cloud file test
                {
                        // 0 bad containers
                        "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram",
                        "gs://hellbender/test/resources/large/human_g1k_v37.20.21.fasta",
                        List.of(Pair.of(CRAMIssue8768Detector.VERBOSE_ARG_NAME, "false")),
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.false.bug8768.cloud.txt",
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.false.bug8768.cloud.tsv",
                },
        };
    }

    @Test(dataProvider = "cramAnalysisTestCases", groups={"cloud"})
    public void testCramAnalysis(
            final String inputPath,
            final String referencePath, // unused for now
            final List<Pair<String, String>> extraArgs,
            final String expectedTextOutputPath,
            final String expectedTSVOutputPath) throws IOException {
        final File outTextFile = createTempFile("testFileDiagnostics", ".txt");
        final File outTSVFile = expectedTSVOutputPath == null ? null : createTempFile("testFileDiagnostics", ".tsv");
        ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addInput(inputPath);
        argBuilder.addOutput(outTextFile);
        argBuilder.addReference(Paths.get(referencePath));
        if (outTSVFile != null) {
            argBuilder.add(CRAMIssue8768Detector.OUTPUT_TSV__ARG_NAME, outTSVFile.getAbsolutePath());
        }
        if (extraArgs != null) {
            extraArgs.forEach(argPair -> argBuilder.add(argPair.getKey(), argPair.getValue()));
        }
        runCommandLine(argBuilder.getArgsList());

        IntegrationTestSpec.assertEqualTextFiles(outTextFile, new File(expectedTextOutputPath));
        if (outTSVFile != null) {
            IntegrationTestSpec.assertEqualTextFiles(outTSVFile, new File(expectedTSVOutputPath));
        }
    }

}
