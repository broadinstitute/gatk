package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class MyCollectDuplexMetricsTest extends CommandLineProgramTest {
    @Test
    public void test(){
        // final String cloud = "gs://fc-secure-429c9379-aa5e-4884-8c35-7a5b947efc37/9bd19b02-cdcd-4454-b332-1277758b8d4b/GenerateDuplexConsensusBams/26ace499-39e5-4985-ae59-25a2ec788c09/call-FGBioGroupReadsByUmi/A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.bam";
        final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
        final String cloud = "gs://broad-dsde-methods-takuto/liquid-biopsy/benchmark-lab-mixture/grouped-by-UMI/A04_denovo_bloodbiopsy_1pct_rep1.fgbio.groupByUmi.bam";
        final String out = "/Users/tsato/workspace/gatk/tmp/my_duplex_metrics.txt";
        final String intervals = "/dsde/working/tsato/consensus/indel/resources/intervals/v1_intersect_denovo.interval_list";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("I", cloud)
                .add("targets", intervals)
                .add("O", out)
                .add(MyCollectDuplexMetrics.READ_NUMBER_LIMIT_NAME, 10);
        runCommandLine(args, MyCollectDuplexMetrics.class.getSimpleName());
    }

}