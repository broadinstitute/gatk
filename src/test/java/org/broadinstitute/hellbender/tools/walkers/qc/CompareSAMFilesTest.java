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
    final String testBam2 = home + "test/Control_KapaLC_1.transcriptome.duplicate_marked_test.sam";

    final String bam1 = home + "Control_KapaLC_1.duplicate_marked.bam";
    final String bam2 = home + "Control_KapaLC_1.transcriptome.duplicate_marked.bam";

    final String gatkHome = "/Users/tsato/workspace/gatk/tmp/";

    final String ribosomeBamHg19 = "";
    final String ribosomeBamHg38 = "";

    @Test
    public void testRibosome(){
        final File out = new File(gatkHome + "ribosome_test.csv");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", ribosomeBamHg19)
                .add("read2", ribosomeBamHg38)
                .add("output", out.getAbsolutePath())
                .add("mode", "RIBOSOME");
        runCommandLine(args, CompareSAMFiles.class.getSimpleName());
        int d = 3;
    }

    @Test
    public void testSmall(){
        final File out = new File(gatkHome + "test.csv");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", testBam1)
                .add("read2", testBam2)
                .add("output", out.getAbsolutePath())
                .add("mode", "DUPLICATE");
        runCommandLine(args, CompareSAMFiles.class.getSimpleName());
        int d = 3;
    }

    @Test
    public void testLarge(){
        final File out = new File(home + "test/test.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", bam1)
                .add("read2", bam2);
        runCommandLine(args, CompareSAMFiles.class.getSimpleName());
        int d = 3;
    }
}