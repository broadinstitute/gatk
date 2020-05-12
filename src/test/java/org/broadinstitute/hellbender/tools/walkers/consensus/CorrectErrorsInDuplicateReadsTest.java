package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class CorrectErrorsInDuplicateReadsTest extends CommandLineProgramTest {
    /**
     * Step through HC code to my heart's desire
     */
    @Test
    public void testHaplotypeCaller(){
        final String bam = "gs://dsde-palantir/SnapshotExperiment2018/HiSeqBams/NA12878_NA12878_IntraRun_1_SM-G947Y_v1.bam";
        final String out = "/Users/tsato/workspace/gatk/tmp/hc.vcf";
        final String intervals = "/Users/tsato/workspace/HipSTR/resources_ts/hg19.hipstr_reference_test.bed";
        final String hg38 = "/seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg38)
                .add("I", bam)
                .add("L", intervals)
                .add("O", out);
        runCommandLine(args, HaplotypeCaller.class.getSimpleName());

        int d = 3;
    }


}