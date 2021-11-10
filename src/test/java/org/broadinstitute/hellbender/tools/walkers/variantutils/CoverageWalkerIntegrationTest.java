package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class CoverageWalkerIntegrationTest extends CommandLineProgramTest {

    private final String contig = "chr16:10000-10000000";

    /** VCF file path (with index in directory) */
    private final String VariantDataPath = ;
    private final File VariantInputData = new File(VariantDataPath);

    /** Associated bam with sample */
    private final String BamPath =
                "gs://dsde-palantir/SnapshotExperiment2018/NovaSeqCRAMs/NA12878_NA12878_O1D1_SM-G947L_v1_NS.cram";

    /** Other bam to compare to */
//    private final String ExtraBamPath =

    /** Output file path */
    private final String outputPath =
            "<my_output_path>/coverage_locuswalker_test.txt";
    private final File outputFile = new File(outputPath);

    /** Hmer interval list file */

    /** Files relating to reference */
    private final File localRefData = new File("<ref file>");

    @Test
    public void testWalker() {

        final List<String> args = Arrays.asList(
//                "-V", VariantDataPath,
                "-I", BamPath,
//                "-eB", ExtraBamPath,
                "-O", outputFile.getPath(),
                "-R", localRefData.getPath(),
                "-L", contig
        );
        runCommandLine(args);
    }
}