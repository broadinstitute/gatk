package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramUnitTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class UMIReadFilterTest extends CommandLineProgramIntegrationTest {
    @Test
    public void test(){
        final String bam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.tiny.bam";
        final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
        // final String umi = "AGG-TGA"; // this is a good one
        final List<String> umis = Arrays.asList("GAC-ATA", "ATT-TTG", "CTG-GCC", "TTG-TAC");
        for (String umi : umis){
            final File out = new File("/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53."+ umi + ".bam");
            final ArgumentsBuilder args = new ArgumentsBuilder()
                    .addArgument("R", hg19)
                    .addArgument("I", bam)
                    .addArgument("read-filter", "UMIReadFilter")
                    .addArgument("umi", umi)
                    .addArgument("O", out.getAbsolutePath());
            runCommandLine(args, PrintReads.class.getSimpleName());
        }

        int d = 3;
    }

    @Test
    public void printConsensusBamWithSpecificUMI(){
        // What did "GAC-ATA" become?
        final String bam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.tiny.bam";
        final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
        final List<String> umis = Arrays.asList("GAC-ATA", "ATT-TTG", "CTG-GCC", "TTG-TAC");
        for (String umi : umis){
            final File out = new File("/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53."+ umi + ".bam");
            final ArgumentsBuilder args = new ArgumentsBuilder()
                    .addArgument("R", hg19)
                    .addArgument("I", bam)
                    .addArgument("read-filter", "UMIReadFilter")
                    .addArgument("umi", umi)
                    .addArgument("O", out.getAbsolutePath());
            runCommandLine(args, PrintReads.class.getSimpleName());
        }

        int d = 3;
    }
}