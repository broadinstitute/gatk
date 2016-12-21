/*
* Copyright 2012-2016 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.rnaseq;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Test all different parameters.
 */
public class ASEReadCounterIntegrationTest extends WalkerTest {

    @Test
    public void testASEReadCounterWithHighMQ() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -o %s -U ALLOW_N_CIGAR_READS", 1,
                Arrays.asList("eec421405a4a570751821d158734020e"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterWithLowMQ() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 1 -o %s -U ALLOW_N_CIGAR_READS", 1,
                Arrays.asList("4f1144be0bb2c4adeae00625afd04cda"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterWithLowMQNoDedup() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 10 -o %s -U ALLOW_N_CIGAR_READS -drf DuplicateRead", 1,
                Arrays.asList("226021673310f28d6520d7f3cfe8cb4b"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterWithHighMQLowBQ() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -o %s -U ALLOW_N_CIGAR_READS", 1,
                Arrays.asList("f86bf14ca3a2cc0114d6a11de0cd9448"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterWithCountOverlaps() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -o %s -U ALLOW_N_CIGAR_READS -overlap COUNT_FRAGMENTS", 1,
                Arrays.asList("bfaeaaa8c000eca82f703225bf2257a1"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterWithCountReads() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -o %s -U ALLOW_N_CIGAR_READS -overlap COUNT_READS", 1,
                Arrays.asList("a9b420fd12a9162b31a48842c4081a1f"));
        executeTest("test high mq with no read passing", spec);
    }


    @Test
    public void testASEReadCounterMinDepth70() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -o %s -U ALLOW_N_CIGAR_READS -minDepth 70", 1,
                Arrays.asList("79b99bc8a1c1c58ac860f79a11f93086"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterFileFormat() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.SYNONYMOUS_CODING.vcf -mmq 60 -mbq 10 -o %s -U ALLOW_N_CIGAR_READS --outputFormat csv", 1,
                Arrays.asList("2c7c531018ab353e6874ee2da7980986"));
        executeTest("test high mq with no read passing", spec);
    }

    @Test
    public void testASEReadCounterNonRef() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ASEReadCounter -R " + b37KGReference + " -I " + privateTestDir + "NA12878.RNAseq.bam -sites "+privateTestDir +"NA12878.chr20_2444518_2637800.RNAseq.NonRef.vcf -U ALLOW_N_CIGAR_READS", 0,
                UserException.class);
        executeTest("test <NON_REF> sites VCF file", spec);
    }










}
