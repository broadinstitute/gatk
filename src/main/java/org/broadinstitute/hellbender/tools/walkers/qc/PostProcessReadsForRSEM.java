package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.http.annotation.Experimental;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Performs post-processing steps to get a bam aligned to a transcriptome ready for RSEM (https://github.com/deweylab/RSEM)
 *
 *
 * Suppose the read name "Q1" aligns to multiple loci in the transcriptome.
 * STAR aligner outputs the reads in the following order:
 *
 * Q1: First-of-pair (Chr 1:1000)
 * Q1: Second-of-pair (Chr 1:2000)
 * Q1: First-of-pair (Chr 20:5000)
 * Q1: Second-of-pair (Chr 20:6000)
 *
 * This is the format required by RSEM. After query-name sorting the reads for duplicate marking,
 * the reads will be ordered as follows;
 *
 * Q1: First-of-pair (Chr 1:1000)
 * Q1: First-of-pair (Chr 20:5000)
 * Q1: Second-of-pair (Chr 1:2000)
 * Q1: Second-of-pair (Chr 20:6000)
 *
 * That is, all the read1's come first, then read2's.
 *
 * This tool reorders such that the alignment pair appears together as in the first example.
 *
 * Caveat: It may be desirable to remove duplicate reads before running RSEM.
 * This tool does not remove duplicate reads; it assumes they have been removed upstream e.g.
 * MarkDuplicates with REMOVE_DUPLICATES=true
 *
 *
 * Usage Example:
 *
 * gatk PostProcessReadsForRSEM \
 * -I input.bam \
 * -O output.bam
 */
@CommandLineProgramProperties(
        summary = "Reorder reads before running RSEM",
        oneLineSummary = "Reorder reads before running RSEM",
        programGroup = ReadDataManipulationProgramGroup.class
)

@BetaFeature
@DocumentedFeature
public class PostProcessReadsForRSEM extends GATKTool {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public GATKPath outSam;

    SAMFileGATKReadWriter writer;

    // Should use CountingFilter or another existing framework for counting and printing filtered reads
    int totalOutputPair = 0;
    int notBothMapped = 0;
    int chimera = 0;
    int unsupportedCigar = 0;

    @Override
    public boolean requiresReads(){
        return true;
    }

    @Override
    public void onTraversalStart(){
        Utils.nonNull(getHeaderForReads(), "Input bam must have a header");
        if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.queryname){
            throw new UserException("Input must be query-name sorted.");
        }

        writer = createSAMWriter(outSam, true);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);
    }

    @Override
    public void traverse() {
        final CountingReadFilter countingReadFilter = makeReadFilter();
        final Iterator<GATKRead> readIterator = getTransformedReadStream(countingReadFilter).iterator();
        ReadPair currentReadPair = null;

        // Initialize the first ReadPair object
        if (readIterator.hasNext()) {
            currentReadPair = new ReadPair(readIterator.next());
            totalOutputPair++;
        }

        while (readIterator.hasNext()){
            final GATKRead read = readIterator.next();
            if (!currentReadPair.getName().equals(read.getName())){
                // We have gathered all the reads with the same query name (primary, secondary, and supplementary alignments)
                // Write the reads to output file, reordering the reads as needed.
                writeReads(currentReadPair);
                currentReadPair = new ReadPair(read);
            } else {
                currentReadPair.add(read);
            }
            progressMeter.update(read);
        }

        if (currentReadPair != null){
            writeReads(currentReadPair);
        }

    }

    /**
     *  For the list of conditions required by RSEM, see: https://github.com/deweylab/RSEM/blob/master/samValidator.cpp
     */
    public boolean passesRSEMFilter(final GATKRead read1, final GATKRead read2) {
        if (read1 == null || read2 == null){
            logger.warn("read1 or read2 is null. This read will not be output. " + read1.getName());
            return false;
        }
        // If either of the pair is unmapped, throw out.
        // With the STAR argument --quantTranscriptomeBan IndelSoftclipSingleend this should not occur,
        // but we check just to be thorough.
        if (read1.isUnmapped() || read2.isUnmapped()){
            notBothMapped++;
            return false;
        }

        // Chimeric reads are not allowed. Chimeric alignments should be just as suspect (or potentially interesting)
        // in the case of transcriptome as in the genome.
        if (!read1.contigsMatch(read2)){
            chimera++;
            return false;
        }

        // Cigar must have a single operator M.
        if (read1.numCigarElements() != 1 || read2.numCigarElements() != 1){
            unsupportedCigar++;
            return false;
        }

        if (read1.getCigarElement(0).getOperator() != CigarOperator.M ||
                read2.getCigarElement(0).getOperator() != CigarOperator.M){
            unsupportedCigar++;
            return false;
        } else {
            return true;
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

        // Write Primary Reads. If either fails, we discard both, and we don't bother with the secondary alignments.
        if (passesRSEMFilter(firstOfPair, secondOfPair)){
            writer.addRead(firstOfPair);
            writer.addRead(secondOfPair);
            totalOutputPair += 1;

            // Now handle secondary alignments.
            final List<Pair<GATKRead, GATKRead>> secondaryAlignmentPairs = groupSecondaryReads(readPair.getSecondaryAlignments());
            for (Pair<GATKRead, GATKRead> mates : secondaryAlignmentPairs){
                // The pair is either both written or both not written.
                if (passesRSEMFilter(mates.getLeft(), mates.getRight())) {
                        writer.addRead(mates.getLeft());
                        writer.addRead(mates.getRight());
                        totalOutputPair += 1;

                }
            }

            // Supplementary reads are not handled i.e. removed
        }
    }

    /**
     *
     * Reorder the secondary alignments as described above.
     * @param secondaryAlignments may be empty.
     * @return a list of pair of matching reads (i.e. read1 with the corresponding read2)
     */
    private List<Pair<GATKRead, GATKRead>> groupSecondaryReads(final List<GATKRead> secondaryAlignments){
        if (secondaryAlignments.isEmpty()){
            return Collections.emptyList();
        }

        final boolean isRead1 = true;
        final Map<Boolean, List<GATKRead>> groupedByRead1 =
                secondaryAlignments.stream().collect(Collectors.groupingBy(r -> r.isFirstOfPair()));
        final List<GATKRead> read1Reads = groupedByRead1.get(isRead1);
        final List<GATKRead> read2Reads = groupedByRead1.get(!isRead1);
        if(read1Reads.size() != read2Reads.size()){
            logger.warn("Num read1s != num read2s among the secondary alignments; " +
                    secondaryAlignments.get(0).getName());
            return Collections.emptyList();
        }


        final List<Pair<GATKRead, GATKRead>> result = new ArrayList<>(read1Reads.size());
        for (GATKRead read1 : read1Reads){
            final List<GATKRead> read2s = read2Reads.stream()
                    .filter(r -> r.contigsMatch(read1)
                            && r.getStart() == read1.getMateStart()
                            && r.getMateStart() == read1.getStart())
                    .collect(Collectors.toList());
            if (read2s.size() == 1){
                result.add(new ImmutablePair<>(read1, read2s.get(0)));
            } else if (read2s.size() > 1){
                logger.warn("Multiple mates found for the secondary alignment " + read1.getName());
            } else {
                logger.warn("Mate not found for the secondary alignment " + read1.getName());
            }
        }
        return result;
    }

    @Override
    public Object onTraversalSuccess(){
        // Write out the last set of reads
        logger.info("Total read pairs output: " + totalOutputPair);
        logger.info("Read pairs filtered due to unmapped mate: " + notBothMapped);
        logger.info("Read pairs filtered due to chimeric alignment: " + chimera);
        logger.info("Read pairs filtered due to a cigar element not supported by RSEM: " + unsupportedCigar);
        return "SUCCESS";
    }

    @Override
    public void closeTool(){
        if (writer != null){
            writer.close();
        }
    }
}
