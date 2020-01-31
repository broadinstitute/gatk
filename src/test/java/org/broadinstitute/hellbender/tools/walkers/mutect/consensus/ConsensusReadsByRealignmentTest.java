package org.broadinstitute.hellbender.tools.walkers.mutect.consensus;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class ConsensusReadsByRealignmentTest extends CommandLineProgramIntegrationTest {
    final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    // This bam is grouped by UMI by fgbio tool
    // final String bam = "/dsde/working/tsato/consensus/bams/C04_denovo_bloodbiopsy_2-5pct_rep1_chr20_umi_grouped_tiny.bam";
    final String bam_ATT_TTG = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.ATT-TTG.bam";

    @Test
    public void test(){
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", bam_ATT_TTG);
        runCommandLine(args);
    }

    @Test
    public void testHaplotypeCallerPloidyOne() throws IOException {
        // final String umi = "ATT-TTG";
        final String umi = "AGG-TGA";
        // final File vcf = createTempFile(umi, "vcf");
        // final File bamout = createTempFile(umi, "bam");
        final String homeDir = "/dsde/working/tsato/consensus/tp53/gatk-debug/";
        final File vcf = new File(homeDir + umi + ".vcf");
        final File bamout = new File(homeDir + umi + ".bam");
        final String bam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53." + umi + ".bam";
        final File graph = new File(homeDir + umi + ".dot");
        bamout.createNewFile();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", bam)
                .addArgument("ploidy", "1")
                .addArgument("L", "17:7578690-7578720") // deli
                .addArgument("force-active")
                .addArgument("debug-assembly")
                .addArgument("graph-output", graph.getAbsolutePath())
                .addArgument("disable-read-filter", "NotDuplicateReadFilter")
                .addArgument("O", vcf.getAbsolutePath())
                .addArgument("bamout", bamout.getAbsolutePath());
        runCommandLine(args, HaplotypeCaller.class.getSimpleName());
        int d = 3;
    }
}