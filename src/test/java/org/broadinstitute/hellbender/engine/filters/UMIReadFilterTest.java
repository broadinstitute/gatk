package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramUnitTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class UMIReadFilterTest extends CommandLineProgramIntegrationTest {
    final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

    // This is a consensus-called bam
    final String C04Chr20Bam = "/dsde/working/tsato/consensus/bams/C04_denovo_bloodbiopsy_2-5pct_rep1.bqsr.bam";
    final String C04Chr20BamPreConcensus = "/dsde/working/tsato/consensus/bams/Jonna_Grimsby_C04_denovo_bloodbiopsy_2-5pct_rep1.bam";

    @Test
    public void test(){
        final String bam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.tiny.bam";
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
    public void printPreConsensusBamWithSpecificUMI(){
        // What did "GAC-ATA" become?
        final String bam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.tiny.bam";
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
    public void findIndels(){
        final File out = new File("/dsde/working/tsato/consensus/indel/C04_denovo_bloodbiopsy_2-5pct_rep1.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", C04Chr20Bam)
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, Mutect2.class.getSimpleName());

        int d = 3;
        // PCR-free data will be a good truth set?
    }

    final File indelVcf = new File(INDEL_DIR + "C04_denovo_bloodbiopsy_2-5pct_rep1_indel.vcf");
    @Test
    public void selectVariants(){
        final File vcf = new File(INDEL_DIR + "C04_denovo_bloodbiopsy_2-5pct_rep1.vcf");
        final File out = indelVcf;

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("V", vcf.getAbsolutePath())
                .addArgument("O", out.getAbsolutePath())
                .addArgument("select-type-to-include", "INDEL");
        runCommandLine(args, SelectVariants.class.getSimpleName());

        int d = 3;
        // PCR-free data will be a good truth set?
    }

    static final String INDEL_DIR = "/dsde/working/tsato/consensus/indel/";
    final File filteredIndelVcf = new File(INDEL_DIR + "C04_denovo_bloodbiopsy_2-5pct_rep1_indel_filtered.vcf");
    @Test
    public void filterMutectCalls(){
        final File mutectStats = new File(INDEL_DIR + "C04_denovo_bloodbiopsy_2-5pct_rep1.vcf.stats");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("V", indelVcf.getAbsolutePath())
                .addArgument("O", filteredIndelVcf.getAbsolutePath())
                .addArgument("stats", mutectStats.getAbsolutePath());
        runCommandLine(args, FilterMutectCalls.class.getSimpleName());

        int d = 3;
    }

    @Test
    public void printPreAndPostConsensusBamWithIndelsC04(){
        // indel positions: 7,579,644,
        // call mutect, find the indel positions, then maybe create an interval file of these spots with
        // also add the "modulo flipping" of UMI for preconsensus only.
        // data is much, much messier than we thought

        final File C04PreConcensusIndelBam = new File("/dsde/working/tsato/consensus/indel/Jonna_Grimsby_C04_denovo_bloodbiopsy_2-5pct_rep1_indels.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", C04Chr20BamPreConcensus)
                .addArgument("read-filter", "RemoveReadsWithoutIndelsReadFilter")
                .addArgument("standardize-umi", "true")
                .addArgument("L", indelVcf.getAbsolutePath())
                .addArgument("O", C04PreConcensusIndelBam.getAbsolutePath());
        runCommandLine(args, PrintReads.class.getSimpleName());

        // Filter the post-consensus bam too
        final File C04ConcensusIndelBam = new File("/dsde/working/tsato/consensus/indel/C04_denovo_bloodbiopsy_2-5pct_rep1.bqsr.indel.bam");
        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", C04Chr20Bam)
                .addArgument("L", indelVcf.getAbsolutePath())
                .addArgument("O", C04ConcensusIndelBam.getAbsolutePath());
        runCommandLine(args2, PrintReads.class.getSimpleName());

        int d = 3;
    }

    @Test
    public void printPreConsensusBamWithIndelsTp53(){
        // indel positions: 7,579,644,
        // call mutect, find the indel positions, then maybe create an interval file of these spots with
        // also add the "modulo flipping" of UMI
        // data is much, much messier than we thought

        final String bam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.tiny.bam";
        final File out = new File("/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.indels_at_17_7578712.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", bam)
                .addArgument("read-filter", "RemoveReadsWithoutIndelsReadFilter")
                .addArgument("L", "17:7578712-7578712")
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, PrintReads.class.getSimpleName());

        // Filter the post-consensus bam too
        final String postConsensusBam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.bqsr.bam";
        final File out2 = new File("/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.bqsr.at_17_7578712.bam");
        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", postConsensusBam)
                .addArgument("L", "17:7578712-7578712")
                .addArgument("O", out2.getAbsolutePath());
        runCommandLine(args2, PrintReads.class.getSimpleName());

        int d = 3;
    }

    @Test
    public void makeSquinsIndelVcf(){
        final String sequinDir = "/dsde/working/tsato/palantir/Analysis/751_Sequins/";
        final String vcf = sequinDir + "sequin_smallvariants_chrQS-germline.vcf";
        final File out = new File(sequinDir + "sequin_smallvariants_chrQS-germline-indel.vcf");
        final String sequinRef = sequinDir + "hg38-chrQ.fasta";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", sequinRef)
                .addArgument("V", vcf)
                .addArgument("select-type-to-include", "INDEL")
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, SelectVariants.class.getSimpleName());

        int d = 3;
    }
}