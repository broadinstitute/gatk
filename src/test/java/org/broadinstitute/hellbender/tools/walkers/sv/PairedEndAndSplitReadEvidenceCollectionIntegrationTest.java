package org.broadinstitute.hellbender.tools.walkers.sv;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;

import java.util.Arrays;

public class PairedEndAndSplitReadEvidenceCollectionIntegrationTest extends GATKBaseTest {

    IntegrationTestSpec spec = new IntegrationTestSpec(
            " -I " + largeFileTestDir + "NA12878.RNAseq.bam -O %s --process-secondary-alignments",
            Arrays.asList(largeFileTestDir + "expected.NA12878.RNAseq.splitNcigarReads.bam"));

}