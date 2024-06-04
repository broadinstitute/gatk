package org.broadinstitute.hellbender;

import org.apache.commons.lang3.tuple.Pair;
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
        return new Object[][]{
                {
                        //this pathname is embedded in the diagnostics output file, so we use a relative pathname
                        // instead of the named constant NA12878_20_21_WGS_cram in order to avoid test failures
                        // caused by the full pathname varying in different environments
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
        };
    }

    @Test(dataProvider = "fileDiagnosticsTestCases")
    public void testFileDiagnostics(
            final String inputPath,
            final List<Pair<String, String>> extraArgs,
            final String expectedOutputPath) throws IOException {
        final File outFile = createTempFile("testFileDiagnostics", ".txt");
        ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addInput(inputPath);
        argBuilder.addOutput(outFile);
        if (extraArgs != null) {
            extraArgs.forEach(argPair -> argBuilder.add(argPair.getKey(), argPair.getValue()));
        }
        runCommandLine(argBuilder.getArgsList());

        IntegrationTestSpec.assertEqualTextFiles(outFile, new File(expectedOutputPath));
    }
}
