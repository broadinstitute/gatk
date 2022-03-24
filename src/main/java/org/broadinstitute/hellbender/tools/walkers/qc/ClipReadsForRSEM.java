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

    @Argument(fullName = "disable-clipping", optional = true)
    public boolean disableClipping = false;

    // For readability
    public boolean enableClipping = !disableClipping;

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
            progressMeter.update(read);

        }
    }

    public boolean passesRSEMFilter(final GATKRead read1, final GATKRead read2) {
        // If either of the pair is unmapped, throw out
        if (read1.getContig() == null || read2.getContig() == null){
            return false;
        }

        // Chimeric reads are not allowed
        if (!read1.getContig().equals(read2.getContig())){
            return false;
        }

        // Cigar must contain at most two elements
        // The only allowed cigar is entirely M, or M followed by S.
        final Cigar cigar1 = read1.getCigar();
        final List<CigarElement> cigarElements1 = cigar1.getCigarElements();
        final Cigar cigar2 = read2.getCigar();
        final List<CigarElement> cigarElements2 = cigar2.getCigarElements();


        if (cigarElements1.size() != cigarElements2.size()){
            // Cigars must have the same length
            return false;
        } else if (cigarElements1.size() == 1){ // And therefore read2 also has size 1
            return cigarElements1.get(0).getOperator() == CigarOperator.M && cigarElements2.get(0).getOperator() == CigarOperator.M;
        } else if (cigarElements1.size() == 2){
            // The only allowed cigars are M followed by S in read1 and S followed by M, and vice versa
            if ((cigarElements1.get(0).getOperator() == CigarOperator.M && cigarElements1.get(1).getOperator() == CigarOperator.S) ||
                    (cigarElements1.get(0).getOperator() == CigarOperator.S && cigarElements1.get(1).getOperator() == CigarOperator.M)){
                // Now check that e.g. 100M46S is paired with 46S100M
                // We don't require the exact match in the sizes of the operators (for now). M
                return cigarElements1.get(0).getOperator() == cigarElements2.get(1).getOperator() &&
                        cigarElements1.get(1).getOperator() == cigarElements2.get(0).getOperator();
            } else {
                return false;
            }


        } else {
            return false;
        }

    }

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
        if (passesRSEMFilter(firstOfPair, secondOfPair)){
            final boolean needsClippingPrimary = enableClipping && needsClipping(firstOfPair, secondOfPair);
            writer.addRead(needsClippingPrimary ? clipRead(firstOfPair) : firstOfPair);
            writer.addRead(needsClippingPrimary ? clipRead(secondOfPair) : secondOfPair);
            outputReadCount += 2;

            final List<Pair<GATKRead, GATKRead>> mateList = groupSecondaryReads(readPair.getSecondaryAlignments());
            for (Pair<GATKRead, GATKRead> mates : mateList){
                // The pair is either both written or both not written
                if (passesRSEMFilter(mates.getLeft(), mates.getRight())){
                    final boolean needsClippingSecondary = enableClipping && needsClipping(mates.getLeft(), mates.getRight());
                    writer.addRead(needsClippingSecondary ? clipRead(mates.getLeft()) : mates.getLeft());
                    writer.addRead(needsClippingSecondary ? clipRead(mates.getRight()) : mates.getRight());
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
            final Optional<GATKRead> read2 = read2Reads.stream().filter(r ->
                            r.getStart() == read1.getMateStart() && r.getMateStart() == read1.getStart()).findFirst();
            if (read2.isPresent()){
                result.add(new ImmutablePair<>(read1, read2.get()));
            } else {
                logger.warn("Mate not found for the secondary alignment " + read1.getName());
            }
        }
        return result;
    }

    // Contract: call passesRSEMFilter on the read pair before calling
    private boolean needsClipping(final GATKRead read1, final GATKRead read2){
        if (read1.getFragmentLength() != -1 * read2.getFragmentLength()){
            logger.warn(read1.getName() + ": Fragment lengths must be negative of each other but got " +
                    read1.getFragmentLength() + ", " + read2.getFragmentLength());
            return false;
        }

        // If only one cigar element, then read1 and read2 are both 146M since we've already run RSEM read filter. No need to clip in this case.
        if (read1.getCigarElements().size() == 1){
            return false;
        }

        return Math.abs(read1.getFragmentLength()) < readLength;
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
        final int length = clippedRead.getLength();

        // For RSEM, remove H from the cigar
        final List<CigarElement> matchCigarElement =  read.getCigarElements().stream().filter(ce -> ce.getOperator() == CigarOperator.M).collect(Collectors.toList());
        Utils.validate(matchCigarElement.size() == 1, "There must be a singl match element but got: " + matchCigarElement);
        // This commented version is the correct way. But sometimes the cigar and read length don't match (a bug in hardClipAdaptorSequence())
        // final CigarElement matchCigarElem = matchCigarElement.get(0);
        final CigarElement matchCigarElem = new CigarElement(clippedRead.getLength(), CigarOperator.M);
        clippedRead.setCigar(new CigarBuilder().add(matchCigarElem).make());
        // This could be off by one, but as long as we get the reads out, we should be ok.
        // Remember this is just a proof of concept; we can fine tune it later, as long as we can verify that we can run RSEM on it.
        Utils.validate(clippedRead.getLength() == matchCigarElem.getLength(),
                "length of cigar operator and read must match but got: " + clippedRead.getLength() + ", " + matchCigarElem.getLength());
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
