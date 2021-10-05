package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.consensus.DownsampleByDuplicateSet;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

public class CompareSAMFilesTest extends CommandLineProgramTest {
    final String home = "/Volumes/dsde_working/tsato/hydro.gen/Analysis/874_twist_RNA/transcriptome/compare_sam_files/";
    final String testBam1 = home + "test/Control_KapaLC_1.duplicate_marked_test.sam";
    final String testBam2 = home + "test/Control_KapaLC_1.transcriptome.dulicate_marked_test.sam";

    final String bam1 = home + "Control_KapaLC_1.duplicate_marked.bam";
    final String bam2 = home + "Control_KapaLC_1.transcriptome.dulicate_marked.bam";


    @Test
    public void test(){
        final File out = new File(home + "test/test.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", bam1)
                .add("read2", bam2);
        runCommandLine(args, CompareSAMFiles.class.getSimpleName());
        int d = 3;
    }
}