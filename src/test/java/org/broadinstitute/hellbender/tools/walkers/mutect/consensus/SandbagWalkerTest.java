package org.broadinstitute.hellbender.tools.walkers.mutect.consensus;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.InferOriginalRead;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

public class SandbagWalkerTest extends CommandLineProgramTest {
    private final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

    @Test
    public void test(){
        final String complexCigar = "/Users/tsato/workspace/gatk/src/test/java/org/broadinstitute/hellbender/tools/walkers/mutect/complex_cigar.bam";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", complexCigar);
        runCommandLine(args, SandbagWalker.class.getSimpleName());
        int d = 3;
    }


}