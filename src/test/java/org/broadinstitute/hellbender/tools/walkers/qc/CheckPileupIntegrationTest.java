/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.qc;

import org.testng.annotations.Test;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;

import java.util.Collections;

/**
 * Run validating pileup across a set of core data as proof of the integrity of the GATK core.
 *
 * Tests both types of old-school pileup formats (basic and consensus).
 *
 * @author mhanna, vdauwera
 * @version 0.2
 */
public class CheckPileupIntegrationTest extends WalkerTest {
    /**
     * This test runs on a consensus pileup containing 10-column lines for SNPs and 13-column lines for indels
     */
    @Test(enabled = true)
    public void testEcoliConsensusPileup() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T CheckPileup" +
                " -I " + validationDataLocation + "MV1994.selected.bam" +
                " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                " --pileup:SAMPileup "+ validationDataLocation + "MV1994.selected.pileup" +
                " -S SILENT -nt 8",0, Collections.<String>emptyList());
        executeTest("testEcoliConsensusPileup",spec);
    }

    /**
     * This test runs on a basic pileup containing 6-column lines for all variants  TODO
     */
    @Test
    public void testEcoliBasicPileup() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T CheckPileup" +
                        " -I " + validationDataLocation + "MV1994.selected.bam" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " --pileup:SAMPileup "+ validationDataLocation + "MV1994.basic.pileup" +
                        " -L Escherichia_coli_K12:1-49" +
                        " -S SILENT -nt 8",0, Collections.<String>emptyList());
        executeTest("testEcoliBasicPileup",spec);
    }
}
