package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.common.annotations.VisibleForTesting;
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
 * Realigns soft-clipped reads. Intended for use with short-read Dragen v3.7.8 alignments.
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

    @Argument(doc="Minimum length of soft clips required for realignment.",
            fullName = "min-clipped-length",
            optional = true,
            minValue = 0)
    public int minSoftClipLength = 0;

    @ArgumentCollection
    public SoftClipRealignmentArgumentCollection softClipRealignmentArgs = new SoftClipRealignmentArgumentCollection();

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
        logger.info("Found " + softclippedReadNames.size() + " soft-clipped reads / read pairs");
        try (final SubsettingRealignmentEngine engine = new SubsettingRealignmentEngine(
                softClipRealignmentArgs.indexImage.toString(), getHeaderForReads(),
                softClipRealignmentArgs.bufferSize, softClipRealignmentArgs.bwaThreads)) {
            logger.info("Subsetting soft-clipped reads and mates...");
            forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                    engine.addRead(read, r -> softclippedReadNames.contains(r.getName()))
            );
            logger.info("Realigning and merging with unclipped reads...");
            final SAMFileGATKReadWriter writer = createSAMWriter(output, true);
            engine.alignAndMerge().forEachRemaining(writer::addRead);
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
    static boolean isValidSoftClip(final GATKRead read, final int minSoftClipLength) {
        return ReadFilterLibrary.PRIMARY_LINE.test(read)
                && read.getCigarElements().stream().anyMatch(c -> c.getOperator() == CigarOperator.SOFT_CLIP
                && c.getLength() >= minSoftClipLength);
    }
}
