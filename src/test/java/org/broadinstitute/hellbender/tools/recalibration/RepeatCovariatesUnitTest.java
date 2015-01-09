/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.RepeatLengthCovariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.RepeatUnitAndLengthCovariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.RepeatUnitCovariate;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class RepeatCovariatesUnitTest {

    RepeatLengthCovariate rlCovariate;
    RepeatUnitCovariate ruCovariate;
    RepeatUnitAndLengthCovariate rurlCovariate;
    RecalibrationArgumentCollection RAC;



    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        rlCovariate = new RepeatLengthCovariate();
        ruCovariate = new RepeatUnitCovariate();
        rurlCovariate = new RepeatUnitAndLengthCovariate();
        rlCovariate.initialize(RAC);
        ruCovariate.initialize(RAC);
        rurlCovariate.initialize(RAC);
    }

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }


    @Test
    public void testFindNumberOfRepetitions() {
        // First, test logic to compute number of repetitions of a substring on a given string.
        int result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACAC".getBytes(), true);
        Assert.assertEquals(2, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACAC".getBytes(), true);
        Assert.assertEquals(4, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACACGT".getBytes(), true);
        Assert.assertEquals(4, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "GTACACACAC".getBytes(), true);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGT".getBytes(), true);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GCAGCAGTAGGGTGTACACACAC".getBytes(), true);
        Assert.assertEquals(1, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GTAGGGTGTACACACACGCAGCAT".getBytes(), true);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGTGTACACACACGCAGCAGCA".getBytes(), true);
        Assert.assertEquals(0, result);
        // Same tests but looking backward on string
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACAC".getBytes(), false);
        Assert.assertEquals(2, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACAC".getBytes(), false);
        Assert.assertEquals(4, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACACGT".getBytes(), false);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "GTACACACAC".getBytes(), false);
        Assert.assertEquals(4, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGT".getBytes(), false);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GCAGCAGTAGGGTGTACACACAC".getBytes(), false);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GTAGGGTGTACACACACGCAGCAT".getBytes(), false);
        Assert.assertEquals(0, result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGTGTACACACACGCAGCAGCA".getBytes(), false);
        Assert.assertEquals(3, result);

        // test logic to get repeat unit and number of repeats from covariate value
        final String[] repUnits = new String[]{"AG","CCG","TCCA","T"};
        for (String ru : repUnits) {
            for (int k=1; k < 10; k++) {
                Pair<String,Integer> pair = RepeatLengthCovariate.getRUandNRfromCovariate(String.format("%s%d", ru, k));
                Assert.assertEquals(pair.getRight().intValue(), k);
                Assert.assertEquals(pair.getLeft(), ru);
            }
        }

    }

    /**
     * Build synthetic reads with random content made up of tandem repeats, record computed Repeat Unit and # repeats and see if
     * they match with read context
     */
    @Test
    public void testManyObservations() {
        final int NUM_UNITS = 10;
        final int MAX_REPEAT_UNIT_LENGTH = RAC.MAX_STR_UNIT_LENGTH;
        final int MAX_NUM_REPETITIONS = RAC.MAX_REPEAT_LENGTH;
        final int NUM_TEST_CASES = 100;

        Random random = new Random();

        for (int r = 0; r < NUM_TEST_CASES; r++) {
            final StringBuilder sb = new StringBuilder();
            // for each unit, generate a repeat unit at random with given random length
            final ArrayList<String> repeatUnits = new ArrayList<String>();
            final ArrayList<Integer> numsRepetitions = new ArrayList<Integer>();
            for (int n=0; n < NUM_UNITS; n++) {
                final int repLength = 1+random.nextInt(MAX_REPEAT_UNIT_LENGTH);
                final String repeatUnit = getRandomBases(repLength);
                final int numRepetitions = 1+random.nextInt(MAX_NUM_REPETITIONS);

                // log for comparison with covariate
                numsRepetitions.add(numRepetitions);
                repeatUnits.add(repeatUnit);

                for (int k=0; k < numRepetitions; k++)
                    sb.append(repeatUnit);

            }

            final String readBases = sb.toString();
            final int readLength = readBases.length();

            final byte[] readQuals = new byte[readLength];
            Arrays.fill(readQuals, (byte) 30);
            final SAMRecord read = ArtificialSAMUtils.createArtificialRead(readBases.getBytes(), readQuals, readLength + "M");

            Covariate[] requestedCovariates = new Covariate[3];
            requestedCovariates[0] = rlCovariate;
            requestedCovariates[1] = ruCovariate;
            requestedCovariates[2] = rurlCovariate;
            ReadCovariates rc = RecalUtils.computeCovariates(read, requestedCovariates);

            // check that the length is correct
            Assert.assertEquals(rc.getMismatchesKeySet().length, readLength);
            Assert.assertEquals(rc.getInsertionsKeySet().length, readLength);
            Assert.assertEquals(rc.getDeletionsKeySet().length, readLength);

            for (int offset = 0; offset < readBases.length(); offset++) { // recalibrate all bases in the read
                // check RepeatLength
                final String rlValM = rlCovariate.formatKey(rc.getMismatchesKeySet(offset)[0]);
                final String rlValI = rlCovariate.formatKey(rc.getInsertionsKeySet(offset)[0]);
                final String rlValD = rlCovariate.formatKey(rc.getDeletionsKeySet(offset)[0]);
                // check RepeatUnit
                final String ruValM = ruCovariate.formatKey(rc.getMismatchesKeySet(offset)[1]);
                final String ruValI = ruCovariate.formatKey(rc.getInsertionsKeySet(offset)[1]);
                final String ruValD = ruCovariate.formatKey(rc.getDeletionsKeySet(offset)[1]);
                // check RepeatUnitAndLength
                final String rurlValM = rurlCovariate.formatKey(rc.getMismatchesKeySet(offset)[2]);
                final String rurlValI = rurlCovariate.formatKey(rc.getInsertionsKeySet(offset)[2]);
                final String rurlValD = rurlCovariate.formatKey(rc.getDeletionsKeySet(offset)[2]);
                // check all 3 values are identical
                Assert.assertEquals(rlValD, rlValI);
                Assert.assertEquals(rlValM, rlValI);
                Assert.assertEquals(ruValD, ruValI);
                Assert.assertEquals(ruValM, ruValI);
                Assert.assertEquals(rurlValD, rurlValI);
                Assert.assertEquals(rurlValM, rurlValI);


                int fw = GATKVariantContextUtils.findNumberOfRepetitions(ruValM.getBytes(), readBases.substring(offset + 1, readLength).getBytes(), true);
                int bw = GATKVariantContextUtils.findNumberOfRepetitions(ruValM.getBytes(), readBases.substring(0, offset + 1).getBytes(), false);
                Assert.assertEquals(Math.min(fw + bw, RAC.MAX_REPEAT_LENGTH), (int) Integer.valueOf(rlValM));
            }

        }
    }

    /**
     * Returns random bases of given length
     * @param length  required length
     * @return        given random string
     */
    String getRandomBases(final int length) {
        byte[] bases = new byte[length];
        Random ran = new Random();
        for (int i=0; i < length; i++ ) {
            int idx = ran.nextInt(4);
            bases[i] = BaseUtils.baseIndexToSimpleBase(idx);
        }
        return new String(bases);
    }


}
