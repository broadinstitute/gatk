package org.broadinstitute.hellbender.tools.evoquer;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

public class EvoquerIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testChr20Dalio3Exomes() throws IOException {
        final String projectID = "broad-dsp-spec-ops";
        final String datasetMapString = "chr20   joint_genotyping_dalio_updated_chr20_3  pet vet";
        final String interval = "chr20:10413398"; // "chr20:10644282" "chr20:10000000-11000000"; //"chr20:157484"; // "chr20:100000-500000";
        final String outputVCF = "/Users/droazen/src/hellbender/evoquer_testChr20Dalio3Exomes.vcf";

        final File datasetMapFile = createTempFile("testChr20Dalio3Exomes", ".dataset_map");
        try ( final PrintWriter writer = new PrintWriter(datasetMapFile) ) {
            writer.println(datasetMapString);
        }

        final String[] args = {
            "--project-id", projectID,
            "--dataset-map", datasetMapFile.getAbsolutePath(),
            "-R", hg38Reference,
            "-L", interval,
            "-O", outputVCF,
            "--run-query-only", "false",
            "--disable-gnarly-genotyper", "false"
        };

        runCommandLine(args);
    }
}
