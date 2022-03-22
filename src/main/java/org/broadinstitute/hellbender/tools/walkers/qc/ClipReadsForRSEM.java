package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

// Should be a GATKTool. ReadWalker filters etc...we don't want any bells and whistles
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReadDataManipulationProgramGroup.class // Sato: Right program group?
)
public class ClipReadsForRSEM extends GATKTool {
    public File outputTable = new File("duplicate_comaprison_counts.csv");

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File outSam;

    @Argument(fullName = "read-length")
    public int readLength = 146;

    PeekableIterator<GATKRead> read1Iterator;
    SAMFileGATKReadWriter writer;

    ReadPair currentReadPair;

    // Check
    int outputReadCount = 0;
    int supplementaryReadCount = 0;
    int totalReadCount = 0;


    @Override
    public void onTraversalStart(){
        read1Iterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());
        // If not query name sorted, then how is this file sorted?
//        if (directlyAccessEngineReadsDataSource().getHeader().getSortOrder() != SAMFileHeader.SortOrder.queryname){
//            throw new UserException("Input must be query-name sorted.");
//        }

        writer = createSAMWriter(new GATKPath(outSam.getAbsolutePath()), true);
        if (!read1Iterator.hasNext()){
            throw new UserException("Input has no reads");
        }

        currentReadPair = new ReadPair(read1Iterator.next());
        outputReadCount++;
    }

    @Override
    public void traverse() {
        while (read1Iterator.hasNext()){
            final GATKRead read = read1Iterator.next();
            if (!currentReadPair.getQueryName().equals(read.getName())){
                // We have gathered all the reads with the same query name (primary, secondary, and supplementary alignments)
                // Write the reads to output file, reordering the reads as needed.
                writeReads(currentReadPair);
                currentReadPair = new ReadPair(read);
            } else {
                currentReadPair.add(read);
            }

        }
    }

    // Wrap this in a read filter
    public boolean passesRSEMFilter(final GATKRead read) {
        // Cigar must contain at most two elements
        // The only allowed cigar is entirely M, or M followed by S.
        final Cigar cigar = read.getCigar();
        final List<CigarElement> cigarElements = cigar.getCigarElements();

        if (cigarElements.size() == 1 && cigarElements.get(0).getOperator() == CigarOperator.M) { // e.g. 146M
            return true;
        } else {
            // e.g. 100M46S or 46S100M. Else false.
            return cigarElements.size() == 2 &&
                    cigarElements.stream().allMatch(ce -> ce.getOperator() == CigarOperator.M ||
                            ce.getOperator() == CigarOperator.S);
        }
    }
    /** Write reads in this order
     *  1. Primary. First of pair, second of pair
     *  2. For each secondary. First of pair, second of pair.
     *  **/
    private void writeReads(final ReadPair readPair){
        // Update read by side effect
        final GATKRead firstOfPair = readPair.getFirstOfPair();
        final GATKRead secondOfPair = readPair.getSecondOfPair();

        // Write Primary Reads. If either fails, we discard both, and we don't bother with the secondary alignments
        // First, check that the read pair passes the RSEM criteria. If not, no need to bother clipping.
        // The only acceptable cigar are 1) 146M, or 2) 142M4S with (insert size) < (read length)
        if (passesRSEMFilter(firstOfPair) && passesRSEMFilter(secondOfPair)){
            writer.addRead(needsClipping(firstOfPair) ? clipRead(firstOfPair) : firstOfPair);
            writer.addRead(needsClipping(secondOfPair) ? clipRead(secondOfPair) : secondOfPair);
            outputReadCount += 2;

            final List<Pair<GATKRead, GATKRead>> mateList = groupSecondaryReads(readPair.getSecondaryAlignments());
            for (Pair<GATKRead, GATKRead> mates : mateList){
                // The pair is either both written or both not written
                if (passesRSEMFilter(mates.getLeft()) && passesRSEMFilter(mates.getRight())){
                    writer.addRead(needsClipping(mates.getLeft()) ? clipRead(mates.getLeft()) : mates.getLeft());
                    writer.addRead(needsClipping(mates.getRight()) ? clipRead(mates.getRight()) : mates.getRight());
                    outputReadCount += 2;
                }
            }

            // Ignore supplementary reads
        }
    }

    /**
     * Contract: secondaryAligments may be empty.
     */
    private List<Pair<GATKRead, GATKRead>> groupSecondaryReads(List<GATKRead> secondaryAlignments){
        if (secondaryAlignments.isEmpty()){
            return Collections.emptyList();
        }

        final boolean READ1 = true;
        final Map<Boolean, List<GATKRead>> groupdbyRead1 =
                secondaryAlignments.stream().collect(Collectors.groupingBy(r -> r.isFirstOfPair()));
        final List<GATKRead> read1Reads = groupdbyRead1.get(READ1);
        final List<GATKRead> read2Reads = groupdbyRead1.get(!READ1);
        Utils.validate(read1Reads.size() == read2Reads.size(), "By assumption we must have the same number of read1s and read2s among the secondary alignments");

        // The pairs is (read1, read2)
        final List<Pair<GATKRead, GATKRead>> result = new ArrayList<>(read1Reads.size());
        for (GATKRead read1 : read1Reads){
            final int mateStart = read1.getMateStart();
            final Optional<GATKRead> read2 = read2Reads.stream().filter(r -> r.getStart() == mateStart).findFirst();
            if (read2.isPresent()){
                result.add(new ImmutablePair<>(read1, read2.get()));
            } else {
                logger.warn("Mate not found for a secondary alignment " + read1.getName());
            }
        }
        return result;
    }

    private boolean needsClipping(final GATKRead read){
        // Also check CIGAR?
        return read.getFragmentLength() < readLength;
    }

    // TODO: replace with ReadClipper.hardClipAdaptorSequence(read)
    private GATKRead clipRead(final GATKRead read){
        // Probably don't need to update mate start (check though). Or start, for that matter
        final Cigar cigar = read.getCigar();
        final byte[] readBases = read.getBases();
        final byte[] quals = read.getBaseQualities();

        final GATKRead clippedRead = ReadClipper.hardClipAdaptorSequence(read);
        final byte[] clippedReadBases = clippedRead.getBases();
        final byte[] clippedQuals = clippedRead.getBaseQualities(); // length = 123...that ok? Should be 124?

        // For RSEM, remove H from the cigar
        final List<CigarElement> matchCigarElement =  read.getCigarElements().stream().filter(ce -> ce.getOperator() == CigarOperator.M).collect(Collectors.toList());
        Utils.validate(matchCigarElement.size() == 1, "There must be a singl match element but got: " + matchCigarElement);
        clippedRead.setCigar(new CigarBuilder().add(matchCigarElement.get(0)).make());

        return clippedRead;
    }



    @Override
    public Object onTraversalSuccess(){
        // Write out the last set of reads
        writeReads(currentReadPair);
        return "SUCCESS";
    }

    @Override
    public void closeTool(){
        writer.close();
    }
}
