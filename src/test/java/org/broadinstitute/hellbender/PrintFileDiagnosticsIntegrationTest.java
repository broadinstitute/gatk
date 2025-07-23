package org.broadinstitute.hellbender;

import htsjdk.beta.plugin.IOUtils;
import htsjdk.io.IOPath;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class PrintFileDiagnosticsIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "fileDiagnosticsTestCases")
    public Object[][] getFileDiagnosticsTestCases() {
        // the pathnames used by the diagnostics tool wind up embedded in the diagnostics output file, so for these
        // tests use just a relative pathname as input (instead of the named constants, i.e., NA12878_20_21_WGS_cram,
        // which are full path names) in order to avoid test failures caused by the full pathname varying in
        // different environments, i.e. in CI
        return new Object[][]{
                {
                        "src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.0.samtools.cram",
                        List.of(Pair.of("count-limit", "10")),
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.txt"
                },
                {
                        "src/test/resources/org/broadinstitute/hellbender/metrics/analysis/MeanQualityByCycle/example_pfFail_reads.bam.bai",
                        null,
                        "src/test/resources/filediagnostics/example_pfFail_reads.bam.bai.txt"
                },
                {
                        "src/test/resources/org/broadinstitute/hellbender/engine/cram_with_crai_index.cram.crai",
                        null,
                        "src/test/resources/filediagnostics/cram_with_crai_index.cram.crai.txt"
                },
                {
                    // cram file that uses all the new 3.1 codecs (fqzcomp, name tok, ransNx16, and adaptive arithmetic)
                        "src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.1.samtools.archive.cram",
                        List.of(Pair.of("count-limit", "20")),
                        "src/test/resources/filediagnostics/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.1.samtools.archive.cram.txt"
                }
        };
    }

    @Test(dataProvider = "fileDiagnosticsTestCases")
    public void testFileDiagnostics(
            final String inputPath,
            final List<Pair<String, String>> extraArgs,
            final String expectedOutputPath) throws IOException {
        final IOPath outFile = IOUtils.createTempPath("testFileDiagnostics", ".txt");
        runFileDiagnosticsTool(new GATKPath(inputPath), extraArgs, outFile);
        IntegrationTestSpec.assertEqualTextFiles(outFile.toPath().toFile(), new File(expectedOutputPath));
    }

    private void runFileDiagnosticsTool(
            final IOPath inputPath,
            final List<Pair<String, String>> extraArgs,
            final IOPath outputPath) {
        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addInput(inputPath.getRawInputString());
        argBuilder.addOutput(outputPath.getRawInputString());
        if (extraArgs != null) {
            extraArgs.forEach(argPair -> argBuilder.add(argPair.getKey(), argPair.getValue()));
        }
        runCommandLine(argBuilder.getArgsList());
    }

}
