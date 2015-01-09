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

import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.ReadClipper;

import java.util.ArrayList;

public class ContextCovariate implements StandardCovariate {
    private final static Logger logger = LogManager.getLogger(ContextCovariate.class);

    private int mismatchesContextSize;
    private int indelsContextSize;

    private int mismatchesKeyMask;
    private int indelsKeyMask;

    private static final int LENGTH_BITS = 4;
    private static final int LENGTH_MASK = 15;

    // the maximum context size (number of bases) permitted; we need to keep the leftmost base free so that values are
    // not negative and we reserve 4 more bits to represent the length of the context; it takes 2 bits to encode one base.
    static final private int MAX_DNA_CONTEXT = 13;
    private byte LOW_QUAL_TAIL;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        indelsContextSize = RAC.INDELS_CONTEXT_SIZE;

        logger.info("\t\tContext sizes: base substitution model " + mismatchesContextSize + ", indel substitution model " + indelsContextSize);

        if (mismatchesContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("mismatches_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, mismatchesContextSize));
        if (indelsContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("indels_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, indelsContextSize));

        LOW_QUAL_TAIL = RAC.LOW_QUAL_TAIL;
        
        if (mismatchesContextSize <= 0 || indelsContextSize <= 0)
            throw new UserException(String.format("Context size must be positive, if you don't want to use the context covariate, just turn it off instead. Mismatches: %d Indels: %d", mismatchesContextSize, indelsContextSize));

        mismatchesKeyMask = createMask(mismatchesContextSize);
        indelsKeyMask = createMask(indelsContextSize);
    }

    @Override
    public void recordValues(final SAMRecord read, final ReadCovariates values) {

        // store the original bases and then write Ns over low quality ones
        final byte[] originalBases = read.getReadBases().clone();
        // Write N's over the low quality tail of the reads to avoid adding them into the context
        final SAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
        
        final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
        byte[] bases = clippedRead.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        final ArrayList<Integer> mismatchKeys = contextWith(bases, mismatchesContextSize, mismatchesKeyMask);
        final ArrayList<Integer> indelKeys = contextWith(bases, indelsContextSize, indelsKeyMask);

        final int readLength = bases.length;

        // this is necessary to ensure that we don't keep historical data in the ReadCovariates values
        // since the context covariate may not span the entire set of values in read covariates
        // due to the clipping of the low quality bases
        if ( readLength != originalBases.length ) {
            // don't both zeroing out if we are going to overwrite the whole array
            for ( int i = 0; i < originalBases.length; i++ )
                // this base has been clipped off, so zero out the covariate values here
                values.addCovariate(0, 0, 0, i);
        }

        for (int i = 0; i < readLength; i++) {
            final int readOffset = (negativeStrand ? readLength - i - 1 : i);
            final int indelKey = indelKeys.get(i);
            values.addCovariate(mismatchKeys.get(i), indelKey, indelKey, readOffset);
        }

        // put the original bases back in
        read.setReadBases(originalBases);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public String formatKey(final int key) {
        if (key == -1) // this can only happen in test routines because we do not propagate null keys to the csv file
            return null;

        return contextFromKey(key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return keyFromContext((String) value);
    }

    private static int createMask(final int contextSize) {
        int mask = 0;
        // create 2*contextSize worth of bits
        for (int i = 0; i < contextSize; i++)
            mask = (mask << 2) | 3;
        // shift 4 bits to mask out the bits used to encode the length
        return mask << LENGTH_BITS;
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *
     * @param bases       the bases in the read to build the context from
     * @param contextSize context size to use building the context
     * @param mask        mask for pulling out just the context bits
     */
    private static ArrayList<Integer> contextWith(final byte[] bases, final int contextSize, final int mask) {

        final int readLength = bases.length;
        final ArrayList<Integer> keys = new ArrayList<Integer>(readLength);

        // the first contextSize-1 bases will not have enough previous context
        for (int i = 1; i < contextSize && i <= readLength; i++)
            keys.add(-1);

        if (readLength < contextSize)
            return keys;

        final int newBaseOffset = 2 * (contextSize - 1) + LENGTH_BITS;

        // get (and add) the key for the context starting at the first base
        int currentKey = keyFromContext(bases, 0, contextSize);
        keys.add(currentKey);

        // if the first key was -1 then there was an N in the context; figure out how many more consecutive contexts it affects
        int currentNPenalty = 0;
        if (currentKey == -1) {
            currentKey = 0;
            currentNPenalty = contextSize - 1;
            int offset = newBaseOffset;
            while (bases[currentNPenalty] != 'N') {
                final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentNPenalty]);
                currentKey |= (baseIndex << offset);
                offset -= 2;
                currentNPenalty--;
            }
        }

        for (int currentIndex = contextSize; currentIndex < readLength; currentIndex++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentIndex]);
            if (baseIndex == -1) { // ignore non-ACGT bases
                currentNPenalty = contextSize;
                currentKey = 0; // reset the key
            } else {
                // push this base's contribution onto the key: shift everything 2 bits, mask out the non-context bits, and add the new base and the length in
                currentKey = (currentKey >> 2) & mask;
                currentKey |= (baseIndex << newBaseOffset);
                currentKey |= contextSize;
            }

            if (currentNPenalty == 0) {
                keys.add(currentKey);
            } else {
                currentNPenalty--;
                keys.add(-1);
            }
        }

        return keys;
    }

    public static int keyFromContext(final String dna) {
        return keyFromContext(dna.getBytes(), 0, dna.length());
    }

    /**
     * Creates a int representation of a given dna string.
     *
     * @param dna    the dna sequence
     * @param start  the start position in the byte array (inclusive)
     * @param end    the end position in the array (exclusive)
     * @return the key representing the dna sequence
     */
    private static int keyFromContext(final byte[] dna, final int start, final int end) {

        int key = end - start;
        int bitOffset = LENGTH_BITS;
        for (int i = start; i < end; i++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
            if (baseIndex == -1) // ignore non-ACGT bases
                return -1;
            key |= (baseIndex << bitOffset);
            bitOffset += 2;
        }
        return key;
    }

    /**
     * Converts a key into the dna string representation.
     *
     * @param key    the key representing the dna sequence
     * @return the dna sequence represented by the key
     */
    public static String contextFromKey(final int key) {
        if (key < 0)
            throw new GATKException("dna conversion cannot handle negative numbers. Possible overflow?");

        final int length = key & LENGTH_MASK; // the first bits represent the length (in bp) of the context
        int mask = 48; // use the mask to pull out bases
        int offset = LENGTH_BITS;

        StringBuilder dna = new StringBuilder();
        for (int i = 0; i < length; i++) {
            final int baseIndex = (key & mask) >> offset;
            dna.append((char)BaseUtils.baseIndexToSimpleBase(baseIndex));
            mask = mask << 2; // move the mask over to the next 2 bits
            offset += 2;
        }

        return dna.toString();
    }

    @Override
    public int maximumKeyValue() {
        // the maximum value is T (11 in binary) for each base in the context
        int length = Math.max(mismatchesContextSize, indelsContextSize);  // the length of the context
        int key = length;
        int bitOffset = LENGTH_BITS;
        for (int i = 0; i <length ; i++) {
            key |= (3 << bitOffset);
            bitOffset += 2;
        }
        return key;
    }
}
