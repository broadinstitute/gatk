package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormalsCreator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Create an allelic panel of of normals (for use in correcting reference bias in the allele-fraction model),
 * given het pulldowns for samples that are part of the panel.  Output to both HDF5 and TSV files is supported.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates an allelic panel of normals, given het pulldowns for samples that are part of the panel, " +
                "and writes it to an HDF5 file (and optionally, a TSV file).  Samples can differ from those " +
                "used to create the coverage PoN but should be representative of the same sequencing process.",
        oneLineSummary = "(Experimental) Creates an allelic panel of normals",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class CreateAllelicPanelOfNormals extends CommandLineProgram {
    protected static final String SITE_FREQUENCY_THRESHOLD_SHORT_NAME = "f";
    protected static final String SITE_FREQUENCY_THRESHOLD_LONG_NAME = "siteFrequencyThreshold";

    protected static final String TSV_OUTPUT_FILE_SHORT_NAME = "tsvOut";
    protected static final String TSV_OUTPUT_FILE_LONG_NAME = "tsvOutput";

    @Argument(
            doc = "Input pulldown files for all samples in the panel of normals.  " +
                    "Duplicate samples are not removed.",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false
    )
    protected List<File> inputFiles = new ArrayList<>();

    @Argument(
            doc = "Output HDF5 file.  The allelic panel of normals will be written to this file " +
                    "(which will be overwritten, if it exists).",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outputFile;

    @Argument(
            doc = "Output TSV file.  The allelic panel of normals will be written as tab-separated values to this file.",
            shortName = TSV_OUTPUT_FILE_SHORT_NAME,
            fullName = TSV_OUTPUT_FILE_LONG_NAME,
            optional = true
    )
    protected File outputTSVFile = null;

    @Argument(
            doc = "Site frequency threshold.  Sites that appear in a fraction of the samples " +
                    "that is strictly less than this threshold are not included in the panel.",
            shortName = SITE_FREQUENCY_THRESHOLD_SHORT_NAME,
            fullName = SITE_FREQUENCY_THRESHOLD_LONG_NAME,
            optional = true
    )
    protected double siteFrequencyThreshold = 0.25;

    @Override
    protected Object doWork() {
        validateArguments();

        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        logger.info("Starting allelic panel of normals creation...");
        final AllelicPanelOfNormals allelicPoN = new AllelicPanelOfNormalsCreator(inputFiles).create(siteFrequencyThreshold);
        logger.info("Allelic panel of normals created.");

        logger.info("Writing allelic panel of normals to output HDF5 file...");
        allelicPoN.write(outputFile, HDF5File.OpenMode.CREATE);
        logger.info("Allelic panel of normals written to " + outputFile + ".");

        if (outputTSVFile != null) {
            allelicPoN.write(outputTSVFile);
            logger.info("Allelic panel of normals written as tab-separated values to " + outputTSVFile + ".");
        }

        return "SUCCESS";
    }

    private void validateArguments() {
        inputFiles.stream().forEach(IOUtils::canReadFile);
        Utils.validateArg(0. < siteFrequencyThreshold && siteFrequencyThreshold <= 1., "Site frequency must be in (0, 1].");
    }
}
