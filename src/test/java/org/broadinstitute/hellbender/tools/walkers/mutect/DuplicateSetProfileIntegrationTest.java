package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Path;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

// (ts 3/16/20): can we delete this class?
public class DuplicateSetProfileIntegrationTest extends CommandLineProgramTest {
    @Test()
    public void test(){
        final String bam = "/dsde/working/tsato/consensus/bams/C04_denovo_bloodbiopsy_2-5pct_rep1_chr20_umi_grouped.bam";

        final List<String> args = Arrays.asList("-I", bam, "-R", b37_reference_20_21);
        runCommandLine(args);
    }


}