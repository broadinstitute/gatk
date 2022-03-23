package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

/**
 * This should eventually be turned into a biread walker, a subclass of WalkerBase.
 * Copied the logic from ComapreSamFiles---refactor.
 *
 * What differentiates this tool from MergeBamAlignment is that
 * - The tags from the unaligned reads are added to the aligned reads, as opposed to the other way.
 *
 *
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "s",
        programGroup = ShortVariantDiscoveryProgramGroup.class // sato: change
)
public class AddUMITagFromUnmappedBam extends GATKTool {
    @Argument(fullName = "unmapped-sam", doc = "unmapped sam file containing the read tag of interest")
    public GATKPath unmappedSamFile;
    ReadsDataSource unmappedSam;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, doc = "")
    public File outSamFile;
    SAMFileGATKReadWriter writer;

    PeekableIterator<GATKRead> alignedSamIterator;
    PeekableIterator<GATKRead> unmappedSamIterator;

    GATKRead currentRead1;
    GATKRead currentRead2;
    ReadQueryNameComparator queryNameComparator;
    QuerynameSetComparison queryNameSetComparison;

    public void onTraversalStart(){
        queryNameComparator = new ReadQueryNameComparator();
        alignedSamIterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());
        SAMFileHeader.SortOrder so1 = directlyAccessEngineReadsDataSource().getHeader().getSortOrder();
        Utils.validate(so1 == SAMFileHeader.SortOrder.queryname, "aligned sam must be sorted by queryname");

        unmappedSam = new ReadsPathDataSource(unmappedSamFile.toPath());
        SAMFileHeader.SortOrder so2 = unmappedSam.getHeader().getSortOrder();
        // Not sure if unmapped bam has sort order field set---check later.
        unmappedSamIterator = new PeekableIterator<>(unmappedSam.iterator());

//        if (alignedSamIterator.hasNext()){
//            currentRead1 = alignedSamIterator.next();
//        } else {
//            throw new UserException("aligned sam iterator is empty.");
//        }
//
        // Initialize read2. Read1 will be initialized in traverse()
        if (unmappedSamIterator.hasNext()){
            currentRead2 = unmappedSamIterator.next();
        } else {
            throw new UserException("unaligned sam iterator is empty.");
        }

        writer = createSAMWriter(new GATKPath(outSamFile.getAbsolutePath()), false);
    }

    final static String UMI_TAG = "RX";

    @Override
    public void traverse() {
        // will need to use this later
        while (alignedSamIterator.hasNext()){
            currentRead1 = alignedSamIterator.next();
            int diff = queryNameComparator.compareReadNames(currentRead1, currentRead2);

            if (diff == 0) {
                // The query names match:
                final String umiTag = currentRead2.getAttributeAsString(UMI_TAG);
                Utils.nonNull(umiTag, "umiTag is empty: read " + currentRead2.getName());
                currentRead1.setAttribute(UMI_TAG, umiTag);
                writer.addRead(currentRead1);

                progressMeter.update(currentRead1);
                continue;
            } else if (diff > 0){
                // Read1 is ahead. Play unmapped reads foward until it catches up.
                while (unmappedSamIterator.hasNext()){
                    currentRead2 = unmappedSamIterator.next();
                    diff = queryNameComparator.compareReadNames(currentRead1, currentRead2);
                    if (diff > 0){
                        continue;
                    } else if (diff == 0){
                        // caught up: star moving the aligned reads forward
                        final String umiTag = currentRead2.getAttributeAsString(UMI_TAG);
                        Utils.nonNull(umiTag, "umiTag is empty: read " + currentRead2.getName());
                        currentRead1.setAttribute(UMI_TAG, umiTag);
                        writer.addRead(currentRead1);
                        break;
                    } else {
                        throw new IllegalStateException("Aligned read is lexicographically smaller than the unmapped read: " +
                                "aligned read = " + currentRead1.getName() + ", unmapped read = " + currentRead2.getName());
                    }
                }
            } else {
                throw new IllegalStateException("Aligned read is lexicographically smaller than the unmapped read: " +
                        "aligned read = " + currentRead1.getName() + ", unmapped read = " + currentRead2.getName());
            }
        }
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

