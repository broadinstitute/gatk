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

package org.broadinstitute.hellbender.tools.recalibration.covariates;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.NGSPlatform;
import org.broadinstitute.hellbender.utils.SequencerFlowClass;

/**
 * The Cycle covariate.
 * For Solexa the cycle is simply the position in the read (counting backwards if it is a negative strand read)
 * For 454 the cycle is the TACG flow cycle, that is, each flow grabs all the TACG's in order in a single cycle
 * For example, for the read: AAACCCCGAAATTTTTACTG
 * the cycle would be 11111111222333333344
 * For SOLiD the cycle is a more complicated mixture of ligation cycle and primer round
 */

public class CycleCovariate implements StandardCovariate {

    private int MAXIMUM_CYCLE_VALUE;
    public static final int CUSHION_FOR_INDELS = 4;
    private String default_platform = null;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        this.MAXIMUM_CYCLE_VALUE = RAC.MAXIMUM_CYCLE_VALUE;

        if (RAC.DEFAULT_PLATFORM != null && !NGSPlatform.isKnown(RAC.DEFAULT_PLATFORM))
            throw new UserException.CommandLineException("The requested default platform (" + RAC.DEFAULT_PLATFORM + ") is not a recognized platform.");

        if (RAC.DEFAULT_PLATFORM != null)
            default_platform = RAC.DEFAULT_PLATFORM;
    }

    // Used to pick out the covariate's value from attributes of the read
    @Override
    public void recordValues(final SAMRecord read, final ReadCovariates values) {
        final int readLength = read.getReadLength();
        final NGSPlatform ngsPlatform = default_platform == null ? NGSPlatform.fromRead(read) : NGSPlatform.fromReadGroupPL(default_platform);

        // Discrete cycle platforms
        if (ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE) {
            final int readOrderFactor = read.getReadPairedFlag() && read.getSecondOfPairFlag() ? -1 : 1;
            final int increment;
            int cycle;
            if (read.getReadNegativeStrandFlag()) {
                cycle = readLength * readOrderFactor;
                increment = -1 * readOrderFactor;
            }
            else {
                cycle = readOrderFactor;
                increment = readOrderFactor;
            }

            final int MAX_CYCLE_FOR_INDELS = readLength - CUSHION_FOR_INDELS - 1;
            for (int i = 0; i < readLength; i++) {
                final int substitutionKey = keyFromCycle(cycle);
                final int indelKey = (i < CUSHION_FOR_INDELS || i > MAX_CYCLE_FOR_INDELS) ? -1 : substitutionKey;
                values.addCovariate(substitutionKey, indelKey, indelKey, i);
                cycle += increment;
            }
        }

        // Flow cycle platforms
        else if (ngsPlatform.getSequencerType() == SequencerFlowClass.FLOW) {

            final byte[] bases = read.getReadBases();

            // Differentiate between first and second of pair.
            // The sequencing machine cycle keeps incrementing for the second read in a pair. So it is possible for a read group
            // to have an error affecting quality at a particular cycle on the first of pair which carries over to the second of pair.
            // Therefore the cycle covariate must differentiate between first and second of pair reads.
            // This effect can not be corrected by pulling out the first of pair and second of pair flags into a separate covariate because
            //   the current sequential model would consider the effects independently instead of jointly.
            final boolean multiplyByNegative1 = read.getReadPairedFlag() && read.getSecondOfPairFlag();

            int cycle = multiplyByNegative1 ? -1 : 1; // todo -- check if this is the right behavior for mate paired reads in flow cycle platforms.

            // BUGBUG: Consider looking at degradation of base quality scores in homopolymer runs to detect when the cycle incremented even though the nucleotide didn't change
            // For example, AAAAAAA was probably read in two flow cycles but here we count it as one
            if (!read.getReadNegativeStrandFlag()) { // Forward direction
                int iii = 0;
                while (iii < readLength) {
                    while (iii < readLength && bases[iii] == (byte) 'T') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    while (iii < readLength && bases[iii] == (byte) 'A') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    while (iii < readLength && bases[iii] == (byte) 'C') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    while (iii < readLength && bases[iii] == (byte) 'G') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    if (iii < readLength) {
                        if (multiplyByNegative1)
                            cycle--;
                        else
                            cycle++;
                    }
                    if (iii < readLength && !BaseUtils.isRegularBase(bases[iii])) {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }

                }
            }
            else { // Negative direction
                int iii = readLength - 1;
                while (iii >= 0) {
                    while (iii >= 0 && bases[iii] == (byte) 'T') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    while (iii >= 0 && bases[iii] == (byte) 'A') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    while (iii >= 0 && bases[iii] == (byte) 'C') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    while (iii >= 0 && bases[iii] == (byte) 'G') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    if (iii >= 0) {
                        if (multiplyByNegative1)
                            cycle--;
                        else
                            cycle++;
                    }
                    if (iii >= 0 && !BaseUtils.isRegularBase(bases[iii])) {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                }
            }
        }

        // Unknown platforms
        else {
            throw new UserException("The platform (" + read.getReadGroup().getPlatform()
                    + ") associated with read group " + read.getReadGroup()
                    + " is not a recognized platform. Allowable options are " + NGSPlatform.knownPlatformsString());
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Integer.parseInt(str);
    }

    @Override
    public String formatKey(final int key) {
        int cycle = key >> 1; // shift so we can remove the "sign" bit
        if ( (key & 1) != 0 ) // is the last bit set?
            cycle *= -1; // then the cycle is negative
        return String.format("%d", cycle);
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? keyFromCycle(Integer.parseInt((String) value)) : keyFromCycle((Integer) value);
    }

    @Override
    public int maximumKeyValue() {
        return (MAXIMUM_CYCLE_VALUE << 1) + 1;
    }

    private int keyFromCycle(final int cycle) {
        // no negative values because values must fit into the first few bits of the long
        int result = Math.abs(cycle);
        if ( result > MAXIMUM_CYCLE_VALUE )
            throw new UserException("The maximum allowed value for the cycle is " + MAXIMUM_CYCLE_VALUE + ", but a larger cycle (" + result + ") was detected.  Please use the --maximum_cycle_value argument to increase this value (at the expense of requiring more memory to run)");

        result = result << 1; // shift so we can add the "sign" bit
        if ( cycle < 0 )
            result++; // negative cycles get the lower-most bit set
        return result;
    }
}