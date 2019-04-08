package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

public class CollectChimeraMapIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test(){
        final String file = "/Users/tsato/workspace/gatk/sequins_only.bam";
        final String sequenceRef = "/humgen/gsa-hpprojects/dev/tsato/palantir/Analysis/751_Sequins/hg38-chrQ.fasta";

        final String[] args = {
                "-R", sequenceRef,
                "-I", file,
                "-O", "counts.tsv"
        };

        runCommandLine(args);

        int d = 3;
    }
}