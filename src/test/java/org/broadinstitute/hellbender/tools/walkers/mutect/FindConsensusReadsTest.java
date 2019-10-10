package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class FindConsensusReadsTest extends CommandLineProgramIntegrationTest {
    @Test
    public void test(){
        final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
        // This bam is grouped by UMI by fgbio tool
        final String bam = "/dsde/working/tsato/consensus/bams/C04_denovo_bloodbiopsy_2-5pct_rep1_chr20_umi_grouped_tiny.bam";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", bam);
        runCommandLine(args);
    }
}