package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * This tool takes a pair of SAM files sharing the same read names (e.g.
 * SAM files before and after alignment) and transfer read tags from one SAM
 * to the other.
 *
 * This situation may happen if for instance an unaligned bam is converted to fastq files,
 * and some read tags get lost during this conversion. The converting to fastq is sometimes unavoidable,
 * due to the fact that some tools, such as the adapter clipping tools, are written for fastq only,
 * while others are written specifically for SAM files.
 *
 * This tools behaves similarly to Picard {@link picard.sam.MergeBamAlignment} (MBA). The difference is that whereas
 * MBA merges the alignment information to the unaligned bam, TransferReadTags uses the aligned bam as the
 * base and adds the read tags from the unaligned bam to the aligned bam.
 *
 * Currently, the tool is implemented for the specific case of transfering UMI read tags (RX)
 * from an unaligned bam.
 */
@CommandLineProgramProperties(
        summary = "Incorporate read tags in a SAM file to that of a matching SAM file",
        oneLineSummary = "Incorporate read tags in a SAM file to that of a matching SAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public class TransferReadTags extends GATKTool {
    @Argument(fullName = "unmapped-sam", doc = "query-name sorted unmapped sam file containing the read tag of interest")
    public GATKPath unmappedSamFile;
    ReadsDataSource unmappedSam;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, doc = "")
    public File outSamFile;
    SAMFileGATKReadWriter writer;

    @Argument(fullName = "read-tags", doc = "read tag names to transfer, must be present in every read in the unmapped sam")
    public List<String> readTags = new ArrayList<>();


    private PeekableIterator<GATKRead> alignedSamIterator;
    private PeekableIterator<GATKRead> unmappedSamIterator;

    private GATKRead currentTargetRead;
    private GATKRead currentUnmappedRead;
    private ReadQueryNameComparator queryNameComparator;

    public void onTraversalStart(){
        Utils.nonEmpty(readTags, "read tags may not be empty");
        queryNameComparator = new ReadQueryNameComparator();
        alignedSamIterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());
        final SAMFileHeader.SortOrder sortOrderAlignedReads = directlyAccessEngineReadsDataSource().getHeader().getSortOrder();
        Utils.validate(sortOrderAlignedReads == SAMFileHeader.SortOrder.queryname, "aligned sam must be sorted by queryname");

        unmappedSam = new ReadsPathDataSource(unmappedSamFile.toPath());

        // We would check that the unmapped bam is query-name sorted, but
        // the SortOrder field is often not populated. Thus we simply assume that
        // the unmapped sam is query-name sorted and rely on the traversal to throw an
        // error when it is not.
        unmappedSamIterator = new PeekableIterator<>(unmappedSam.iterator());

        // Initialize the current unmapped read. The counterpart for the aligned read will be initialized in traverse()
        if (unmappedSamIterator.hasNext()){
            currentUnmappedRead = unmappedSamIterator.next();
        } else {
            throw new UserException("unmapped sam iterator is empty.");
        }

        writer = createSAMWriter(new GATKPath(outSamFile.getAbsolutePath()), false);
    }

    /**
     * The traversal assumes that the aligned (target) reads is a subset of
     * the unmapped reads. The tool exists when there exists a read in the
     * aligned read that is not present in the unmapped read.
     */
    @Override
    public void traverse() {
        // will need to use this later
        while (alignedSamIterator.hasNext()){
            currentTargetRead = alignedSamIterator.next();
            int diff = queryNameComparator.compareReadNames(currentTargetRead, currentUnmappedRead);

            if (diff == 0) {
                // The query names match.
                final GATKRead updatedRead = updateReadTags(currentTargetRead, currentUnmappedRead);
                writer.addRead(updatedRead);
                progressMeter.update(currentTargetRead);
                continue;
            } else if (diff > 0){
                // target read is ahead; play unmapped reads forward until it catches up.
                while (unmappedSamIterator.hasNext()){
                    currentUnmappedRead = unmappedSamIterator.next();
                    diff = queryNameComparator.compareReadNames(currentTargetRead, currentUnmappedRead);
                    if (diff > 0){
                        continue;
                    } else if (diff == 0){
                        // caught up: start moving the aligned reads forward
                        final GATKRead updatedRead = updateReadTags(currentTargetRead, currentUnmappedRead);
                        writer.addRead(updatedRead);
                        break;
                    } else {
                        throw new IllegalStateException("Aligned read is lexicographically smaller than the unmapped read. This tool assumes " +
                                "reads in both input files are query-name sorted lexicographically (i.e. by Picard SortSam but not by samtools sort): " +
                                "aligned read = " + currentTargetRead.getName() + ", unmapped read = " + currentUnmappedRead.getName() + "");
                    }
                }
            } else {
                throw new IllegalStateException("Aligned read is lexicographically smaller than the unmapped read. This tool assumes " +
                        "reads in both input files are query-name sorted lexicographically (i.e. by Picard SortSam but not by samtools sort): " +
                        "aligned read = " + currentTargetRead.getName() + ", unmapped read = " + currentUnmappedRead.getName() + "");
            }
        }
    }

    /**
     *
     * @param targetRead We use this read as the template
     * @param originRead The read from which we extract the requested fields
     * @return A new instance of GATKRead derived from targetRead, updated with the requested fields from originRead
     */
    private GATKRead updateReadTags(final GATKRead targetRead, final GATKRead originRead){
        final GATKRead updatedRead = targetRead.copy();
        for (String tagName : readTags){
            final String tagValue = originRead.getAttributeAsString(tagName);
            Utils.nonNull(tagValue, "The attribute is empty: read " + currentUnmappedRead.getName());
            updatedRead.setAttribute(tagName, tagValue);
        }
        return updatedRead;
    }

    @Override
    public Object onTraversalSuccess(){
        Utils.validate(!alignedSamIterator.hasNext(), "aligned sam iterator has to have iterated through");
        if (writer != null) {
            writer.close();
        }
        return "SUCCESS";
    }
}

