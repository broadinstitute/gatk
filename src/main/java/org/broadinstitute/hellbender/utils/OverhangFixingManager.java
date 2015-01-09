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

package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

import java.util.*;

/**
 * The class manages reads and splices and tries to apply overhang clipping when appropriate.
 * Important note: although for efficiency the manager does try to send reads to the underlying writer in coordinate
 * sorted order, it does NOT guarantee that it will do so in every case!  So unless there's a good reason not to,
 * methods that instantiate this manager should pass in a writer that does not assume the reads are pre-sorted.
 */
public class OverhangFixingManager {

    protected static final Logger logger = LogManager.getLogger(OverhangFixingManager.class);
    private static final boolean DEBUG = false;

    // how many reads should we store in memory before flushing the queue?
    private final int MAX_RECORDS_IN_MEMORY;

    // how many mismatches do we tolerate in the overhangs?
    private final int MAX_MISMATCHES_IN_OVERHANG;

    // how many bases do we tolerate in the overhang before deciding not to clip?
    private final int MAX_BASES_IN_OVERHANG;

    // should we not bother fixing overhangs?
    private final boolean doNotFixOverhangs;

    // where we ultimately write out our records
    private final SAMFileWriter writer;

    // fasta reference reader to check overhanging edges in the exome reference sequence
    private final IndexedFastaSequenceFile referenceReader;

    // the genome loc parser
    private final GenomeLocParser genomeLocParser;

    // the read cache
    private final static int initialCapacity = 5000;
    private PriorityQueue<SplitRead> waitingReads = new PriorityQueue<>(initialCapacity, new SplitReadComparator());

    // the set of current splices to use
    private final Set<Splice> splices = new TreeSet<>(new SpliceComparator());

    protected static final int MAX_SPLICES_TO_KEEP = 1000;


    /**
     *
     * @param writer                   actual writer
     * @param genomeLocParser          the GenomeLocParser object
     * @param referenceReader          the reference reader
     * @param maxRecordsInMemory       max records to keep in memory
     * @param maxMismatchesInOverhangs max number of mismatches permitted in the overhangs before requiring clipping
     * @param maxBasesInOverhangs      max number of bases permitted in the overhangs before deciding not to clip
     * @param doNotFixOverhangs        if true, don't clip overhangs at all
     */
    public OverhangFixingManager(final SAMFileWriter writer,
                                 final GenomeLocParser genomeLocParser,
                                 final IndexedFastaSequenceFile referenceReader,
                                 final int maxRecordsInMemory,
                                 final int maxMismatchesInOverhangs,
                                 final int maxBasesInOverhangs,
                                 final boolean doNotFixOverhangs) {
        this.writer = writer;
        this.genomeLocParser = genomeLocParser;
        this.referenceReader = referenceReader;
        this.MAX_RECORDS_IN_MEMORY = maxRecordsInMemory;
        this.MAX_MISMATCHES_IN_OVERHANG = maxMismatchesInOverhangs;
        this.MAX_BASES_IN_OVERHANG = maxBasesInOverhangs;
        this.doNotFixOverhangs = doNotFixOverhangs;
    }

    public final int getNReadsInQueue() { return waitingReads.size(); }

    /**
     * For testing purposes only
     *
     * @return the list of reads currently in the queue
     */
    public List<SplitRead> getReadsInQueueForTesting() {
        return new ArrayList<>(waitingReads);
    }

    /**
     * For testing purposes only
     *
     * @return the list of splices currently in the queue
     */
    public List<Splice> getSplicesForTesting() {
        return new ArrayList<>(splices);
    }

    /**
     * Add a new observed split to the list to use
     *
     * @param contig  the contig
     * @param start   the start of the split, inclusive
     * @param end     the end of the split, inclusive
     */
    public void addSplicePosition(final String contig, final int start, final int end) {
        if ( doNotFixOverhangs )
            return;

        // is this a new splice?  if not, we are done
        final Splice splice = new Splice(contig, start, end);
        if ( splices.contains(splice) )
            return;

        // initialize it with the reference context
        // we don't want to do this until we know for sure that it's a new splice position
        splice.initialize(referenceReader);

        // clear the set of old split positions seen if we hit a new contig
        final boolean sameContig = splices.isEmpty() || splices.iterator().next().loc.getContig().equals(contig);
        if ( !sameContig )
            splices.clear();

        // run this position against the existing reads
        for ( final SplitRead read : waitingReads )
            fixSplit(read, splice);

        splices.add(splice);

        if ( splices.size() > MAX_SPLICES_TO_KEEP )
            cleanSplices();
    }

    /**
     * Add a read to the manager
     *
     * @param read  the read to add
     */
    public void addRead(final SAMRecord read) {
        if ( read == null ) throw new IllegalArgumentException("read added to manager is null, which is not allowed");

        // if the new read is on a different contig or we have too many reads, then we need to flush the queue and clear the map
        final boolean tooManyReads = getNReadsInQueue() >= MAX_RECORDS_IN_MEMORY;
        final boolean encounteredNewContig = getNReadsInQueue() > 0 && !waitingReads.peek().read.getReferenceIndex().equals(read.getReferenceIndex());

        if ( tooManyReads || encounteredNewContig ) {
            if ( DEBUG ) logger.warn("Flushing queue on " + (tooManyReads ? "too many reads" : ("move to new contig: " + read.getReferenceName() + " from " + waitingReads.peek().read.getReferenceName())) + " at " + read.getAlignmentStart());

            final int targetQueueSize = encounteredNewContig ? 0 : MAX_RECORDS_IN_MEMORY / 2;

            // write the required number of waiting reads to disk
            while ( getNReadsInQueue() > targetQueueSize )
                writer.addAlignment(waitingReads.poll().read);
        }

        final SplitRead splitRead = new SplitRead(read);

        // fix overhangs, as needed
        for ( final Splice splice : splices)
            fixSplit(splitRead, splice);

        // add the new read to the queue
        waitingReads.add(splitRead);
    }

    /**
     * Clean up the list of splices
     */
    private void cleanSplices() {
        final int targetQueueSize = splices.size() / 2;
        final Iterator<Splice> iter = splices.iterator();
        for ( int i = 0; i < targetQueueSize; i++ ) {
            iter.next();
            iter.remove();
        }
    }

    /**
     * Try to fix the given read using the given split
     *
     * @param read        the read to fix
     * @param splice      the split (bad region to clip out)
     */
    private void fixSplit(final SplitRead read, final Splice splice) {
        // if the read doesn't even overlap the split position then we can just exit
        if ( read.loc == null || !splice.loc.overlapsP(read.loc) )
            return;

        if ( isLeftOverhang(read.loc, splice.loc) ) {
            final int overhang = splice.loc.getStop() - read.loc.getStart() + 1;
            if ( overhangingBasesMismatch(read.read.getReadBases(), 0, splice.reference, splice.reference.length - overhang, overhang) ) {
                final SAMRecord clippedRead = ReadClipper.hardClipByReadCoordinates(read.read, 0, overhang - 1);
                read.setRead(clippedRead);
            }
        }
        else if ( isRightOverhang(read.loc, splice.loc) ) {
            final int overhang = read.loc.getStop() - splice.loc.getStart() + 1;
            if ( overhangingBasesMismatch(read.read.getReadBases(), read.read.getReadLength() - overhang, splice.reference, 0, overhang) ) {
                final SAMRecord clippedRead = ReadClipper.hardClipByReadCoordinates(read.read, read.read.getReadLength() - overhang, read.read.getReadLength() - 1);
                read.setRead(clippedRead);
            }
        }
    }

    /**
     * Is this a proper overhang on the left side of the read?
     *
     * @param readLoc    the read's loc
     * @param spliceLoc   the split's loc
     * @return true if it's a left side overhang
     */
    protected static boolean isLeftOverhang(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        return readLoc.getStart() <= spliceLoc.getStop() && readLoc.getStart() > spliceLoc.getStart() && readLoc.getStop() > spliceLoc.getStop();
    }

    /**
     * Is this a proper overhang on the right side of the read?
     *
     * @param readLoc    the read's loc
     * @param spliceLoc   the split's loc
     * @return true if it's a right side overhang
     */
    protected static boolean isRightOverhang(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        return readLoc.getStop() >= spliceLoc.getStart() && readLoc.getStop() < spliceLoc.getStop() && readLoc.getStart() < spliceLoc.getStart();
    }

    /**
     * Are there too many mismatches to the reference among the overhanging bases?
     *
     * @param read                  the read bases
     * @param readStartIndex        where to start on the read
     * @param reference             the reference bases
     * @param referenceStartIndex   where to start on the reference
     * @param spanToTest            how many bases to test
     * @return true if too many overhanging bases mismatch, false otherwise
     */
    protected boolean overhangingBasesMismatch(final byte[] read,
                                               final int readStartIndex,
                                               final byte[] reference,
                                               final int referenceStartIndex,
                                               final int spanToTest) {
        // don't process too small a span, too large a span, or a span that is most of a read
        if ( spanToTest < 1 || spanToTest > MAX_BASES_IN_OVERHANG || spanToTest > read.length / 2 )
            return false;

        int numMismatchesSeen = 0;
        for ( int i = 0; i < spanToTest; i++ ) {
            if ( read[readStartIndex + i] != reference[referenceStartIndex + i] ) {
                if ( ++numMismatchesSeen > MAX_MISMATCHES_IN_OVERHANG )
                    return true;
            }
        }

        // we can still mismatch overall if at least half of the bases mismatch
        return numMismatchesSeen >= ((spanToTest+1)/2);
    }

    /**
     * Close out the manager stream by clearing the read cache
     */
    public void close() {
        // write out all of the remaining reads
        while ( ! waitingReads.isEmpty() )
            writer.addAlignment(waitingReads.poll().read);
    }

    // class to represent the reads with their soft-clip-included GenomeLocs
    public final class SplitRead {

        public SAMRecord read;
        public GenomeLoc loc;

        public SplitRead(final SAMRecord read) {
            setRead(read);
        }

        public void setRead(final SAMRecord read) {
            boolean readIsEmpty = ReadUtils.isEmpty(read);
            if ( !readIsEmpty ) {
                this.read = read;
                if ( ! read.getReadUnmappedFlag() )
                    loc = genomeLocParser.createGenomeLoc(read.getReferenceName(), ReadUtils.getSoftStart(read), ReadUtils.getSoftEnd(read));
            }
        }
    }

    // class to represent the comparator for the split reads
    private final class SplitReadComparator implements Comparator<SplitRead> {

        private final SAMRecordCoordinateComparator readComparator;

        public SplitReadComparator() {
            readComparator = new SAMRecordCoordinateComparator();
        }

        public int compare(final SplitRead read1, final SplitRead read2) {
            return readComparator.compare(read1.read, read2.read);
        }
    }

    // class to represent the split positions
    protected final class Splice {

        public final GenomeLoc loc;
        public byte[] reference;

        public Splice(final String contig, final int start, final int end) {
            loc = genomeLocParser.createGenomeLoc(contig, start, end);
        }

        public void initialize(final IndexedFastaSequenceFile referenceReader) {
            reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        }

        @Override
        public boolean equals(final Object other) {
            return other != null && (other instanceof Splice) && this.loc.equals(((Splice)other).loc);
        }

        @Override
        public int hashCode() {
            return loc.hashCode();
        }
    }

    // class to represent the comparator for the split reads
    private final class SpliceComparator implements Comparator<Splice> {

        public int compare(final Splice position1, final Splice position2) {
            return position1.loc.compareTo(position2.loc);
        }
    }
}
