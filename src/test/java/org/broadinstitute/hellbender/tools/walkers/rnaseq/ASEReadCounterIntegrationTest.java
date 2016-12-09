package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class ASEReadCounterIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testASEReadCounterWithHighMQ() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -O %s -mmq 60 ",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.WithHighMQ.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithLowMQ() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 1 -O %s ",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.WithLowMQ.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithLowMQNoDedup() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 10 -O %s -DF NotDuplicateReadFilter",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.WithLowMQNoDedup.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithHighMQLowBQ() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s ",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.WithHighMQLowBQ.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithCountOverlaps() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s -overlap COUNT_FRAGMENTS",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.WithCountOverlaps.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterWithCountReads() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s -overlap COUNT_READS",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.WithCountReads.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterMinDepth70() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s -minDepth 70",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.MinDepth70.table"));
        spec.executeTest("test high mq with no read passing", this);
    }

    @Test
    public void testASEReadCounterFileFormat() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-R " + b37_reference_20_21 + " -I " + largeFileTestDir + "NA12878.RNAseq.bam -V " + publicTestDir + "NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -O %s --outputFormat CSV",
                Arrays.asList(publicTestDir + "expected.ASEReadCount.FileFormat.table"));
        spec.executeTest("test high mq with no read passing", this);
    }
}
