package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiplePassReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Realigns soft-clipped reads using BWA.
 *
 * <p>
 *
 * </p>
 *
 * <p>
 * The reference is strictly required when handling CRAM files.
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 * Print reads that pass the WellformedReadFilter and that map to chromosome 20.
 * <pre>
 * gatk PrintReads \
 *   -I input.bam \
 *   -L 20 \
 *   -O chr20.bam
 * </pre>
 * Print reads that pass the WellformedReadFilter and the NotDuplicateReadFilter.
 * <pre>
 * gatk PrintReads \
 *   -I input.bam \
 *   --read-filter NotDuplicateReadFilter \
 *   -O filtered.bam
 * </pre>
 *
 * Print reads that pass the WellformedReadFilter and that map to chromosome 20 from a BAM in a google cloud bucket to another cloud bucket.
 * <pre>
 * gatk PrintReads \
 *   -I gs://cloud-bucket/input.bam \
 *   -L 20 \
 *   -O gs://my-gcs-bucket/chr20.bam
 * </pre>
 * Print reads that pass the WellformedReadFilter and that map to chromosome 20 from a BAM in a google cloud bucket to the local system.
 * <pre>
 * gatk PrintReads \
 *   -I gs://cloud-bucket/input.bam \
 *   -L 20 \
 *   -O chr20.bam
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Realigns soft-clipped reads to a given reference.",
        oneLineSummary = "Realigns soft-clipped reads to a given reference.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public final class RealignSoftClippedReads extends MultiplePassReadWalker {

    @Argument(doc="Output bam file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public GATKPath output;

    @Argument(doc="Minimum length of soft clips required for realignment. Setting to 0 will realign all reads.",
            fullName = "min-clipped-length",
            optional = true,
            minValue = 0)
    public int minSoftClipLength = 1;

    @ArgumentCollection
    public SubsettingRealignmentArgumentCollection args = new SubsettingRealignmentArgumentCollection();

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public void traverseReads() {
        logger.info("Identifying soft-clipped reads...");
        final Set<String> softclippedReadNames = new HashSet<>();
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                checkIfClipped(read, softclippedReadNames, minSoftClipLength)
        );
        logger.info("Found " + softclippedReadNames.size() + " soft-clipped reads / read pairs.");
        try (final SubsettingRealignmentEngine engine = new SubsettingRealignmentEngine(args.indexImage.toString(),
                getHeaderForReads(), args.bufferSize, args.bwaThreads, args.keepDuplicateFlag);
             final SAMFileGATKReadWriter writer = createSAMWriter(output, true)) {
            logger.info("Subsetting soft-clipped reads and mates...");
            forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) -> {
                // Clear realignment tags in case the input was previously realigned
                read.clearAttribute(SubsettingRealignmentEngine.REALIGNED_READ_TAG);
                // Submit the read to the engine, selecting soft-clipped reads and their mates
                engine.addRead(read, r -> softclippedReadNames.contains(r.getName()));
            });
            logger.info("Found " + engine.getSelectedReadsCount() + " soft-clipped reads and mates.");
            logger.info("Found " + engine.getNonselectedReadsCount() + " non-clipped reads.");
            logger.info("Realigning and merging reads...");
            engine.alignAndMerge().forEachRemaining(writer::addRead);
            logger.info("Realigned " + engine.getPairedAlignmentReadsCount() + " paired and " +
                    engine.getUnpairedAlignmentReadsCount() + " unpaired reads.");
        }
    }

    /**
     * Gather names of soft-clipped reads
     */
    @VisibleForTesting
    static void checkIfClipped(final GATKRead read, final Set<String> softclippedReadNames, final int minSoftClipLength) {
        if (isValidSoftClip(read, minSoftClipLength)) {
            softclippedReadNames.add(read.getName());
        }
    }

    /**
     * Returns whether the read is non-secondary/non-supplementary and sufficiently soft-clipped
     */
    @VisibleForTesting
    static boolean isValidSoftClip(final GATKRead read, final int minSoftClipLength) {
        return ReadFilterLibrary.PRIMARY_LINE.test(read) && hasMinSoftClip(read, minSoftClipLength);
    }

    private static boolean hasMinSoftClip(final GATKRead read, final int minSoftClipLength) {
        if (minSoftClipLength == 0) {
            return true;
        }
        for (final CigarElement e : read.getCigarElements()) {
            if (e.getOperator() == CigarOperator.SOFT_CLIP && e.getLength() >= minSoftClipLength) {
                return true;
            }
        }
        return false;
    }
}
