package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

public class AddUMITagFromUnmappedBamIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test(){
        final String home = "/Volumes/dsde_working/tsato/hydro.gen/Analysis/874_twist_RNA/add_umi/";

//        final File alignedBam = new File(home + "SM-KYN26_SSIV_SM-LQZZ2_transcriptome.grouped.queryname_sorted.bam");
//        final File unmappedSam = new File(home + "SM-KYN26_SSIV_SM-LQZZ2_UMI_extracted.bam");
//        final File out = new File(home + "SM-KYN26_SSIV_SM-LQZZ2_transcriptome_with_UMI.bam");

        final File alignedBam = new File(home + "SM-LVFDV_15_Min_Low_High_transcriptome.grouped.queryname_sorted.bam");
        final File unmappedSam = new File(home + "SM-LVFDV_15_Min_Low_High_UMI_extracted.bam");
        final File out = new File(home + "SM-LVFDV_15_Min_Low_High_transcriptome_with_UMI.bam");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", alignedBam)
                .add("unmapped-sam", unmappedSam)
                .add("O", out.getAbsolutePath());
        runCommandLine(args, AddUMITagFromUnmappedBam.class.getSimpleName());
        int d = 3;
    }


}