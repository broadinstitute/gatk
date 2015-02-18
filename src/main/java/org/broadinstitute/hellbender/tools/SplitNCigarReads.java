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

package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.bqsr.ReadTransformer;
import org.broadinstitute.hellbender.transformers.NDNCigarReadTransformer;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.OverhangFixingManager;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.sam.CigarUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.*;

/**
 *
 * Splits reads that contain Ns in their cigar string (e.g. spanning splicing events).
 *
 * Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements).
 * The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N
 * (including the Ns) is hard clipped and so on for the rest of the new reads.
 */
@CommandLineProgramProperties(
        usage = "Splits reads that contain Ns in their cigar string (e.g. spanning splicing events).",
        usageShort = "Split Reads with N in Cigar",
        programGroup = ReadProgramGroup.class
)
public class SplitNCigarReads extends CommandLineProgram {

    @Argument(fullName = INPUT_LONG_NAME, shortName= INPUT_SHORT_NAME, doc="The SAM/BAM/CRAM file to read from.")
    public File INPUT;

    @Argument(fullName = OUTPUT_LONG_NAME, shortName = OUTPUT_SHORT_NAME, doc="Write output to this BAM filename instead of STDOUT")
    protected File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.",
            common = true, optional = true)
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    /**
     * This flag tells GATK to refactor cigar string with NDN elements to one element. It intended primarily for use in
     * a RNAseq pipeline since the problem might come up when using RNAseq aligner such as Tophat2 with provided transcriptoms.
     * You should only use this if you know that your reads have that problem.
     */
    @Argument(fullName = "refactor_NDN_cigar_string", shortName = "fixNDN", doc = "refactor cigar string with NDN elements to one element", optional = true)
    public boolean REFACTOR_NDN_CIGAR_READS = false;

    /**
     * For expert users only!  To minimize memory consumption you can lower this number, but then the tool may skip
     * overhang fixing in regions with too much coverage.  Just make sure to give Java enough memory!  4Gb should be
     * enough with the default value.
     */
    @Argument(fullName="maxReadsInMemory", shortName="maxInMemory", doc="max reads allowed to be kept in memory at a time by the BAM writer", optional=true)
    protected int MAX_RECORDS_IN_MEMORY = 150000;

    /**
     * If there are more than this many mismatches within the overhang regions, the whole overhang will get hard-clipped out.
     * It is still possible in some cases that the overhang could get clipped if the number of mismatches do not exceed this
     * value, e.g. if most of the overhang mismatches.
     */
    @Argument(fullName="maxMismatchesInOverhang", shortName="maxMismatches", doc="max number of mismatches allowed in the overhang", optional=true)
    protected int MAX_MISMATCHES_IN_OVERHANG = 1;

    /**
     * If there are more than this many bases in the overhang, we won't try to hard-clip them out
     */
    @Argument(fullName="maxBasesInOverhang", shortName="maxOverhang", doc="max number of bases allowed in the overhang", optional=true)
    protected int MAX_BASES_TO_CLIP = 40;

    @Argument(fullName="doNotFixOverhangs", shortName="doNotFixOverhangs", doc="do not have the walker hard-clip overhanging sections of the reads", optional=true)
    protected boolean doNotFixOverhangs = false;

    /**
     * This stores all of the already-split reads and manages any processing (e.g. clipping overhangs) that happens to them.
     * It will emit reads to the underlying writer as needed so we don't need to worry about any of that in this class.
     */
    protected OverhangFixingManager overhangManager;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileWriter outputWriter = initialize(in);

        final ReadTransformer rnaReadTransform = REFACTOR_NDN_CIGAR_READS ? new NDNCigarReadTransformer() : ReadTransformer.identity();

        StreamSupport.stream(in.spliterator(), false)
                .map(rnaReadTransform)
                .forEach(read -> splitNCigarRead(read, overhangManager));
        overhangManager.close();
        CloserUtil.close(in);
        CloserUtil.close(outputWriter);
        return null;
    }

    private SAMFileWriter initialize(final SamReader in) {
        final SAMFileHeader outputHeader = in.getFileHeader().clone();
        final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, REFERENCE_SEQUENCE);

        try {
            final IndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(REFERENCE_SEQUENCE);
            GenomeLocParser genomeLocParser= new GenomeLocParser(referenceReader.getSequenceDictionary());
            overhangManager = new OverhangFixingManager(outputWriter, genomeLocParser, referenceReader, MAX_RECORDS_IN_MEMORY, MAX_MISMATCHES_IN_OVERHANG, MAX_BASES_TO_CLIP, doNotFixOverhangs);
            return outputWriter;
        } catch (FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(REFERENCE_SEQUENCE, ex);
        }
    }

    /**
     * Goes through the cigar string of the read and create new reads for each consecutive non-N elements (while hard clipping the rest of the read).
     * For example: for a read with cigar '1H2M2D1M2N1M2I1N1M2S' 3 new reads will be created with cigar strings: 1H2M2D1M, 1M2I and 1M2S
     *
     * @param read     the read to split
     */
    public static SAMRecord splitNCigarRead(final SAMRecord read, OverhangFixingManager manager) {
        final int numCigarElements = read.getCigar().numCigarElements();

        int firstCigarIndex = 0;
        for ( int i = 0; i < numCigarElements; i++ ) {
            final CigarElement cigarElement = read.getCigar().getCigarElement(i);
            if (cigarElement.getOperator() == CigarOperator.N) {
                manager.addRead(splitReadBasedOnCigar(read, firstCigarIndex, i, manager));
                firstCigarIndex = i+1;
            }
        }

        // if there are no N's in the read
        if (firstCigarIndex == 0) {
            manager.addRead(read);
        }
        //add the last section of the read: from the last N to the the end of the read
        // (it will be done for all the usual cigar string that does not end with N)
        else if (firstCigarIndex < numCigarElements) {
            manager.addRead(splitReadBasedOnCigar(read, firstCigarIndex, numCigarElements, null));
        }
        return read;
    }

    /**
     * Pull out an individual split position for a read
     *
     * @param read               the read being split
     * @param cigarStartIndex    the index of the first cigar element to keep
     * @param cigarEndIndex      the index of the last cigar element to keep
     * @param forSplitPositions  the manager for keeping track of split positions; can be null
     * @return a non-null read representing the section of the original read being split out
     */
    private static SAMRecord splitReadBasedOnCigar(final SAMRecord read, final int cigarStartIndex, final int cigarEndIndex, final OverhangFixingManager forSplitPositions) {
        int cigarFirstIndex = cigarStartIndex;
        int cigarSecondIndex = cigarEndIndex;

        //in case a section of the read ends or starts with D (for example the first section in 1M1D1N1M is 1M1D), we should trim this cigar element
        // it can be 'if', but it was kept as 'while' to make sure the code can work with Cigar strings that were not "cleaned"
        while(read.getCigar().getCigarElement(cigarFirstIndex).getOperator().equals(CigarOperator.D))
            cigarFirstIndex++;
        while(read.getCigar().getCigarElement(cigarSecondIndex-1).getOperator().equals(CigarOperator.D))
            cigarSecondIndex--;
        if(cigarFirstIndex > cigarSecondIndex)
            throw new IllegalArgumentException("Cannot split this read (might be an empty section between Ns, for example 1N1D1N): "+read.getCigarString());

        // we keep only the section of the read that is aligned to the reference between startRefIndex and stopRefIndex (inclusive).
        // the other sections of the read are clipped:
        final int startRefIndex = read.getUnclippedStart() + CigarUtils.countRefBasesBasedOnCigar(read, 0, cigarFirstIndex); //goes through the prefix of the cigar (up to cigarStartIndex) and move the reference index.
        final int stopRefIndex = startRefIndex + CigarUtils.countRefBasesBasedOnCigar(read,cigarFirstIndex,cigarSecondIndex)-1; //goes through a consecutive non-N section of the cigar (up to cigarEndIndex) and move the reference index.

        if ( forSplitPositions != null ) {
            final String contig = read.getReferenceName();
            final int splitStart = startRefIndex + CigarUtils.countRefBasesBasedOnCigar(read,cigarFirstIndex,cigarEndIndex);  //we use cigarEndIndex instead of cigarSecondIndex so we won't take into account the D's at the end.
            final int splitEnd = splitStart + read.getCigar().getCigarElement(cigarEndIndex).getLength() - 1;
            forSplitPositions.addSplicePosition(contig, splitStart, splitEnd);
        }

        return ReadClipper.hardClipToRegionIncludingClippedBases(read, startRefIndex, stopRefIndex);
    }

}
