package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class VariantQualityWalkerIntegrationTest extends CommandLineProgramTest {
    /** Includes index file in local directory */
    private final File localVariantInputData = new File("/Users/rmagner/Documents/Data/JB_Samples/001029_B0245_NA12878_JBX-045_BC13_E02.g.vcf.gz");

    /** Output file path */
    private final File outputFile = new File( "/Users/rmagner/Documents/Data/JB_Samples/output.tsv");

    /** Files relating to reference */
    private final File localRefData = new File("/Users/rmagner/Documents/Data/Reference/Homo_sapiens_assembly38.fasta");
    private final File hmerRefTable = new File("/Users/rmagner/Documents/Data/Reference/HmerData/hmer_hg38_db.table");

    @Test
    public void testWalker() {
//        final File outputFile = createTempFile("output", ".txt");

        final List<String> args = Arrays.asList(
                "-R", localRefData.getPath(),
                "-O", outputFile.getPath(),
                "-av", localVariantInputData.getPath(),
                "-F", hmerRefTable.getPath(),
                "-L", "chr21");
        runCommandLine(args);
    }
}