package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class PhantomBaseLocusWalkerIntegrationTest extends CommandLineProgramTest {

    private final String contig = "chr18:108074-120000";

    /** VCF file path (with index in directory) */
    private final String VariantDataPath = "gs://fc-secure-b116e047-dc4d-49b4-9661-92c38055e426/66c13f31-1b1f-46b5-9521-31c619111f39/JukeboxVC/526242af-f44b-4426-93dc-29c3112bd81b/call-FilterSymbolicAlleles/001029_B0245_NA12878_JBX-045_BC13_E02.fix.vcf.gz";
    private final File VariantInputData = new File(VariantDataPath);

    /** Associated bam with sample */
    private final String BamPath =
            "gs://fc-secure-0c5a3c69-1298-4e9d-8b33-055d9291ad01/e4009d5c-667b-4be4-af5f-be8a5d74e979/JukeboxVC/8e54e7ee-0099-4609-aa62-55e786e3c998/call-ConvertToCram/002404_NA12878-No_DTT-X0007.cram";
    //        "gs://fc-secure-b116e047-dc4d-49b4-9661-92c38055e426/66c13f31-1b1f-46b5-9521-31c619111f39/JukeboxVC/f36a010b-afba-4fcd-af4b-89f2f2e7b8e9/call-ConvertToCram/001029_B0245_NA12878_JBX-045_BC18_B03.cram";
    //"gs://dsde-palantir/SnapshotExperiment2018/NovaSeqCRAMs/NA12878_NA12878_O1D1_SM-G947L_v1_NS.cram";
//    private final File BamFile = new File(BamPath);

    /** Output file path */
    private final String outputPath =
            "/Users/rmagner/Documents/pb_locuswalker_test.txt";
    //     "/Users/rmagner/Documents/Analysis/PhantomBases/Data/Jukebox_001029_B0245_NA12878_JBX-045_BC18_B03-gt4-" + contig + ".tsv";
    //"/Users/rmagner/Documents/Analysis/PhantomBases/Data/Illumina_NA12878_goldstandard_O1D1_SM-G947L_v1_NS-gt4-" + contig + ".tsv" ;
    private final File outputFile = new File(outputPath);

    /** Hmer interval list file */
    private final File hmerIntervalList = new File("/Users/rmagner/Documents/Data/Reference/HmerData/hmer_hg38.list");

    /** Files relating to reference */
    private final File localRefData = new File("/Users/rmagner/Documents/Data/Reference/Homo_sapiens_assembly38.fasta");
    private final File hmerRefTable = new File("/Users/rmagner/Documents/Data/Reference/HmerData/hmer_hg38_gt4.table");

    @Test
    public void testWalker() {

        final List<String> args = Arrays.asList(
//                "-V", VariantDataPath,
                "-I", BamPath,
                "-O", outputFile.getPath(),
                "-R", localRefData.getPath(),
                "-F", hmerRefTable.getPath(),
                "-L", contig
        );
//                "-L", hmerIntervalList.getPath() );
        runCommandLine(args);
    }
}