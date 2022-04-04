package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Performs post-processing steps to get a bam aligned to a transcriptome ready for RSEM (https://github.com/deweylab/RSEM)
 *
 * ### Task 1. Reordering Reads ###
 *
 * Suppose the queryname "Q1" aligns to multiple loci in the transcriptome.
 * STAR aligner outputs the reads in the following order:
 *
 * Q1: Read1 (Chr 1:1000)
 * Q1: Read2 (Chr 1:2000)
 * Q1: Read1 (Chr 20:5000)
 * Q1: Read2 (Chr 20:6000)
 *
 * This is the format required by RSEM. After query-name sorting the reads for duplicate marking,
 * the reads will be ordered as follows;
 *
 * Q1: Read1 (Chr 1:1000)
 * Q1: Read1 (Chr 20:5000)
 * Q1: Read2 (Chr 1:2000)
 * Q1: Read2 (Chr 20:6000)
 *
 * That is, all the read1's come first, then read2's.
 *
 * This tool reorders such that the alignment pair appears together as in the first example.
 *
 * ### Task 2. Removing Reads ###
 *
 * If requested, this tool also removes duplicate marked reads and MT reads, which can skew gene expression counts.
 *
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReadDataManipulationProgramGroup.class // Sato: Change to QC when the other PR is merged.
)
public class PostProcessReadsForRSEM extends GATKTool {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File outSam;

    @Argument(fullName = "keep-MT-reads")
    public boolean keepMTReads = false;

    @Argument(fullName = "keep-duplicates")
    public boolean keepDuplicates = false;

    PeekableIterator<GATKRead> read1Iterator;
    SAMFileGATKReadWriter writer;

    ReadPair currentReadPair;

    // Use CountingFilter or another existing framework for counting and printing filtered reads
    int totalOutputPair = 0;
    int notBothMapped = 0;
    int chimera = 0;
    int unsupportedCigar = 0;
    int duplicates = 0;
    int mitochondria = 0;

    @Override
    public void onTraversalStart(){
        read1Iterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());

        if (directlyAccessEngineReadsDataSource().getHeader().getSortOrder() != SAMFileHeader.SortOrder.queryname){
            throw new UserException("Input must be query-name sorted.");
        }

        writer = createSAMWriter(new GATKPath(outSam.getAbsolutePath()), true);

        if (!read1Iterator.hasNext()){
            throw new UserException("Input has no reads");
        }

        currentReadPair = new ReadPair(read1Iterator.next());
        totalOutputPair++;
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
        // If either of the pair is unmapped, throw out.
        // With the STAR argument --quantTranscriptomeBan IndelSoftclipSingleend this should not occur,
        // but we check just to be thorough.
        if (read1.getContig() == null || read2.getContig() == null){
            notBothMapped++;
            return false;
        }

        // Chimeric reads are not allowed
        if (!read1.getContig().equals(read2.getContig())){
            chimera++;
            return false;
        }

        // Cigar must have a single operator M.
        final List<CigarElement> cigarElements1 = read1.getCigar().getCigarElements();
        final List<CigarElement> cigarElements2 = read2.getCigar().getCigarElements();

        if (cigarElements1.size() != 1 || cigarElements2.size() != 1){
            unsupportedCigar++;
            return false;
        }

        if (cigarElements1.get(0).getOperator() != CigarOperator.M ||
                cigarElements2.get(0).getOperator() != CigarOperator.M){
            unsupportedCigar++;
            return false;
        } else {
            return true;
        }
    }

    private final static HashSet<String> MT_TRANSCRIPT_IDs = new HashSet<>(Arrays.asList(
            "ENST00000387314.1", "ENST00000389680.2", "ENST00000387342.1", "ENST00000387347.2", "ENST00000386347.1",
            "ENST00000361390.2", "ENST00000387365.1", "ENST00000387372.1", "ENST00000387377.1", "ENST00000361453.3",
            "ENST00000387382.1", "ENST00000387392.1", "ENST00000387400.1", "ENST00000387405.1", "ENST00000387409.1",
            "ENST00000361624.2", "ENST00000387416.2", "ENST00000387419.1", "ENST00000361739.1", "ENST00000387421.1",
            "ENST00000361851.1", "ENST00000361899.2", "ENST00000362079.2", "ENST00000387429.1", "ENST00000361227.2",
            "ENST00000387439.1", "ENST00000361335.1", "ENST00000361381.2", "ENST00000387441.1", "ENST00000387449.1",
            "ENST00000387456.1", "ENST00000361567.2", "ENST00000361681.2", "ENST00000387459.1", "ENST00000361789.2",
            "ENST00000387460.2", "ENST00000387461.2"
    ));
    private boolean mitochondrialRead(final GATKRead read){
        return MT_TRANSCRIPT_IDs.contains(read.getContig());
    }

    private boolean passesAdditionalFilters(final GATKRead r1, final GATKRead r2){
        if (! keepDuplicates){
            // If the read pair is duplicate, then return false
            if (r1.isDuplicate()){
                duplicates++;
                return false;
            }
        }

        if (!keepMTReads){
            // If the read pair is duplicate, then return false
            if (mitochondrialRead(r1) || mitochondrialRead(r2)){
                mitochondria++;
                return false;
            }
        }

        return true;
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
        if (passesRSEMFilter(firstOfPair, secondOfPair) && passesAdditionalFilters(firstOfPair, secondOfPair)){
            writer.addRead(firstOfPair);
            writer.addRead(secondOfPair);
            totalOutputPair += 1;

            // Now handle secondary alignments.
            final List<Pair<GATKRead, GATKRead>> secondaryAlignmentPairs = groupSecondaryReads(readPair.getSecondaryAlignments());
            for (Pair<GATKRead, GATKRead> mates : secondaryAlignmentPairs){
                // The pair is either both written or both not written.
                if (passesRSEMFilter(mates.getLeft(), mates.getRight()) &&
                        passesAdditionalFilters(mates.getLeft(), mates.getRight())) {
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
        final Map<Boolean, List<GATKRead>> groupdByRead1 =
                secondaryAlignments.stream().collect(Collectors.groupingBy(r -> r.isFirstOfPair()));
        final List<GATKRead> read1Reads = groupdByRead1.get(isRead1);
        final List<GATKRead> read2Reads = groupdByRead1.get(!isRead1);
        if(read1Reads.size() != read2Reads.size()){
            logger.warn("Num read1s != num read2s among the secondary alignments; " +
                    secondaryAlignments.get(0).getName());
            return Collections.emptyList();
        }


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

    @Override
    public Object onTraversalSuccess(){
        // Write out the last set of reads
        writeReads(currentReadPair);
        logger.info("Total read pairs output: " + totalOutputPair);
        logger.info("Read pairs filtered due to unmapped mate: " + notBothMapped);
        logger.info("Read pairs filtered due to chimeric alignment: " + chimera);
        logger.info("Read pairs filtered due to a cigar element not supported by RSEM: " + unsupportedCigar);
        return "SUCCESS";
    }

    @Override
    public void closeTool(){
        writer.close();
    }
}
