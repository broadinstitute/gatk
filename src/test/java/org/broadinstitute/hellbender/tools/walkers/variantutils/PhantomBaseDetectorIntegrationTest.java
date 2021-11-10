package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class PhantomBaseDetectorIntegrationTest extends CommandLineProgramTest {

    /** VCF file path (with index in directory) */
    private final String VariantDataPath = "gs://fc-secure-b116e047-dc4d-49b4-9661-92c38055e426/66c13f31-1b1f-46b5-9521-31c619111f39/JukeboxVC/526242af-f44b-4426-93dc-29c3112bd81b/call-FilterSymbolicAlleles/001029_B0245_NA12878_JBX-045_BC13_E02.fix.vcf.gz";
    private final File VariantInputData = new File(VariantDataPath);

    /** Associated bam with sample */
    private final String BamPath = "gs://fc-secure-b116e047-dc4d-49b4-9661-92c38055e426/66c13f31-1b1f-46b5-9521-31c619111f39/JukeboxVC/f36a010b-afba-4fcd-af4b-89f2f2e7b8e9/call-ConvertToCram/001029_B0245_NA12878_JBX-045_BC18_B03.cram";
//    private final File BamFile = new File(BamPath);

    /** Output file path */
    private final String outputPath = "/Users/rmagner/Documents/Data/Jukebox/test.txt" ;
    private final File outputFile = new File(outputPath);

    /** Hmer interval list file */
    private final File hmerIntervalList = new File("/Users/rmagner/Documents/Data/Reference/HmerData/hmer_hg38.list");

    /** Files relating to reference */
    private final File localRefData = new File("/Users/rmagner/Documents/Data/Reference/Homo_sapiens_assembly38.fasta");
    private final File hmerRefTable = new File("/Users/rmagner/Documents/Data/Reference/HmerData/hmer_hg38_gt6.table");

    @Test
    public void testWalker() {

        final List<String> args = Arrays.asList(
//                "-V", VariantDataPath,
                "-I", BamPath,
                "-O", outputFile.getPath(),
                "-R", localRefData.getPath(),
                "-F", hmerRefTable.getPath(),
                "-L", "chr21"
        );
//                "-L", hmerIntervalList.getPath() );
        runCommandLine(args);
    }
}