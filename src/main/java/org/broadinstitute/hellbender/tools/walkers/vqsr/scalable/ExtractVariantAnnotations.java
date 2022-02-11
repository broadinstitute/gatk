package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.List;

/**
 * TODO
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class ExtractVariantAnnotations extends VariantLabeledAnnotationsWalker {

    @Override
    public boolean isExtractVariantsNotOverlappingResources() {
        return false;
    }

    @Override
    public Object onTraversalSuccess() {

        if (dataBatch.size() == 0) {
            throw new GATKException("None of the specified input variants were present in the resource VCFs.");
        }

        for (final String resourceLabel : dataBatch.sortedLabels) {
            logger.info(String.format("Extracted annotations for %d variants labeled as %s.",
                    (int) dataBatch.getData().stream().flatMap(List::stream).mapToDouble(vd -> vd.labels.contains(resourceLabel) ? 1 : 0).sum(),
                    resourceLabel));
        }
        logger.info(String.format("Extracted annotations for %s total variants.", dataBatch.size()));

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }
}