package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

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
public final class ExtractVariantAnnotations extends VariantAnnotationWalker {

    @Override
    public void beforeOnTraversalStart() {
        isExtractAll = false;
    }

    @Override
    public void afterTraversalSuccess() {

        if (data.getData().isEmpty()) {
            throw new GATKException("None of the specified input variants were present in the resource VCFs.");
        }

        // TODO iterate through labels and report totals
//        logger.info(String.format("Extracted annotations for %s truth variants.", dataManager.getData().stream().filter(v -> v.atTruthSite).count()));
        logger.info(String.format("Extracted annotations for %s total variants.", data.getData().size()));

        logger.info("Writing annotations...");
        writeAnnotationsHDF5();

        logger.info("Writing VCF...");
        writeVCF(true, true, false);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));
    }
}