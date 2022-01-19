package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;

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

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output annotations HDF5 file.")
    private File outputAnnotationsHDF5File;

    @Override
    public void beforeOnTraversalStart() {
        isExtractTrainingAndTruthOnly = true;

        // fail early if we cannot create the output file
        if ((outputAnnotationsHDF5File.exists() && !outputAnnotationsHDF5File.canWrite()) ||
                (!outputAnnotationsHDF5File.exists() && !outputAnnotationsHDF5File.getAbsoluteFile().getParentFile().canWrite())) {
            throw new UserException(String.format("Cannot create output file at %s.", outputAnnotationsHDF5File));
        }
    }

    @Override
    public void afterTraversalSuccess() {

        if (dataManager.getData().isEmpty()) {
            throw new GATKException("None of the specified input variants were present in the resource VCFs.");
        }
        logger.info(String.format("Extracted annotations for %s truth variants.", dataManager.getData().stream().filter(v -> v.atTruthSite).count()));
        logger.info(String.format("Extracted annotations for %s total variants.", dataManager.getData().size()));

        logger.info("Writing annotations...");
        writeAnnotationsHDF5(outputAnnotationsHDF5File);
        logger.info(String.format("Annotations written to %s.", outputAnnotationsHDF5File.getAbsolutePath()));

        logger.info(String.format("%s complete.", getClass().getSimpleName()));
    }
}