package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 *
 * Tests all possible (and valid) cigar strings that might contain any cigar elements. It uses a code that were written to test the ReadClipper walker.
 * For valid cigar sting in length 8 there are few thousands options, with N in every possible option and with more than one N (for example 1M1N1M1N1M1N2M).
 * The cigarElements array is used to provide all the possible cigar element that might be included.
 */
public final class SplitNCigarReadsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testSplitsWithOverhangs()  throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -O %s ",
                Arrays.asList(largeFileTestDir + "expected.NA12878.RNAseq.splitNcigarReads.bam"));    //results created using gatk3.46
        spec.executeTest("test splits with overhangs", this);
    }

    @Test
    public void testSplitsWithOverhangsNotClipping() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "--doNotFixOverhangs -R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -O %s ",
                Arrays.asList(largeFileTestDir + "expected.NA12878.RNAseq.splitNcigarReads.doNotFixOverhangs.bam"));   //results created using gatk3.46
        spec.executeTest("test splits with overhangs not clipping", this);
    }

    @Test
    public void testSplitsWithOverhangs0Mismatches() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "--maxMismatchesInOverhang 0 -R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -O %s ",
                Arrays.asList(largeFileTestDir + "expected.NA12878.RNAseq.splitNcigarReads.maxMismatchesInOverhang0.bam"));   //results created using gatk3.46
        spec.executeTest("test splits with overhangs 0 mismatches", this);
    }

    @Test
    public void testSplitsWithOverhangs5BasesInOverhang()  throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "--maxBasesInOverhang 5 -R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -O %s ",
                Arrays.asList(largeFileTestDir + "expected.NA12878.RNAseq.splitNcigarReads.maxBasesInOverhang5.bam"));    //results created using gatk3.46
        spec.executeTest("test splits with overhangs 5 bases in overhang", this);
    }

    @Test
    public void testSplitsFixNDN() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + getTestDataDir() +"/" + "splitNCigarReadsSnippet.bam -O %s -fixNDN",
                Arrays.asList(getTestDataDir() +"/" + "expected.splitNCigarReadsSnippet.splitNcigarReads.fixNDN.bam"));
        spec.executeTest("test fix NDN", this);
    }
}
