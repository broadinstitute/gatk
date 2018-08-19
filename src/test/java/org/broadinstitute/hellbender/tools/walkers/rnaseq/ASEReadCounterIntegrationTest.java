package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class ASEReadCounterIntegrationTest extends CommandLineProgramTest {
    public static final String aseTestDir = toolsTestDir + "walkers/rnaseq/ASEReadCounter/";

    @Test
    public void testASEReadCounterWithHighMQ() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -O %s -mmq 60 ",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.WithHighMQ.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithLowMQ() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 1 -O %s ",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.WithLowMQ.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithLowMQNoDedup() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 10 -O %s -DF NotDuplicateReadFilter",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.WithLowMQNoDedup.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithHighMQLowBQ() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s ",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.WithHighMQLowBQ.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithCountOverlaps() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s -overlap COUNT_FRAGMENTS",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.WithCountOverlaps.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithCountReads() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s -overlap COUNT_READS",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.WithCountReads.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterMinDepth70() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s -min-depth 70",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.MinDepth70.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterFileFormat() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s --output-format CSV",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.FileFormat.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test(expectedExceptions = UserException.class)
    public void testASEReadCounterMultipleContexts() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-R");
        args.add(b37_reference_20_21);
        args.add("-I");
        args.add(largeFileTestDir + "NA12878.RNAseq.bam");
        args.add("-V");
        args.add(aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.MultiContext.vcf");
        args.add("-O");
        args.add(GATKBaseTest.createTempFile("testMultipleContexts", ".csv"));

        runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.class)
    public void testASEReadCounterNonRefAllele() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-R");
        args.add(b37_reference_20_21);
        args.add("-I");
        args.add(largeFileTestDir + "NA12878.RNAseq.bam");
        args.add("-V");
        args.add(aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.NON_REF.vcf");
        args.add("-O");
        args.add(GATKBaseTest.createTempFile("testMultipleContexts", ".csv"));

        runCommandLine(args);
    }

    @Test
    public void testASEReadCounterWarnings() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-R");
        args.add(b37_reference_20_21);
        args.add("-I");
        args.add(largeFileTestDir + "NA12878.RNAseq.bam");
        args.add("-V");
        args.add(aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.warnings.vcf");
        args.add("-O");
        args.add(GATKBaseTest.createTempFile("testMultipleContexts", ".csv"));

        runCommandLine(args);
    }

    @Test
    public void testASEReadCounterImproperPairs() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -V " + aseTestDir + "NA12878.chr20_2444518_2637800.RNAseq.IMPROPER_PAIR.vcf -mmq 60 -mbq 10 -O %s --output-format CSV",
                Arrays.asList(aseTestDir + "expected.ASEReadCount.ImproperPair.table"));
        spec.executeTest("test improper pairs", this);
    }
}
