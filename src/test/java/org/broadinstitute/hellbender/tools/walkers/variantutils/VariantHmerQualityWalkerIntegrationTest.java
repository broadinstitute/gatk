package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class VariantHmerQualityWalkerIntegrationTest extends CommandLineProgramTest {

    private final int contig = 22;
    private final String sequencer = "Jukebox";
    // Illumina sample: "SnapshotExperiment2015_CEUTrio_gvcfs_NA12878.d691bf66-37af-4375-9c09-6a6607f322e8.g.vcf.gz"
    // Jukebox sample: "001029_B0245_NA12878_JBX-045_BC13_E02.g.vcf.gz"
    private final String inputGVCFName = "001029_B0245_NA12878_JBX-045_BC13_E02.g.vcf.gz";

    /** Includes index file in local directory */
    private final String localVariantDataPath = "/Users/rmagner/Documents/Data/" + sequencer + "/" + inputGVCFName;
    private final File localVariantInputData = new File(localVariantDataPath);

    /** Output file path */
//    private final String outputPath = "/Users/rmagner/Documents/Data/" + sequencer + "/VHQW/VHQW_output_chr" + contig + ".tsv";
//    private final String outputPath = "/Users/rmagner/Documents/Data/" + sequencer + "/VHQW/indelsAndSNPs-" + contig + ".tsv";
    private final String outputPath = "/Users/rmagner/Documents/test.tsv";
    private final File outputFile = new File(outputPath);

    /** Files relating to reference */
    private final File localRefData = new File("/Users/rmagner/Documents/Data/Reference/Homo_sapiens_assembly38.fasta");
    private final File hmerRefTable = new File("/Users/rmagner/Documents/Data/Reference/HmerData/hmer_hg38_db.table");

    @Test
    public void testWalker() {
//        final File outputFile = createTempFile("output", ".txt");

        final List<String> args = Arrays.asList(
                "-R", localRefData.getPath(),
                "-O", outputFile.getPath(),
                "-V", localVariantInputData.getPath(),
                "-F", hmerRefTable.getPath(),
                "-L", "chr" + contig );
        runCommandLine(args);
    }
}