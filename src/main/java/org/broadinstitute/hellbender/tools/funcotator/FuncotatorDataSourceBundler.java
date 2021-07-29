package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.codec.digest.MessageDigestAlgorithms;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeter;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeterResults;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * {@link FuncotatorDataSourceBundler} is a tool to download data sources for a specified organism for <b><i>{@link Funcotator}</i></b>.
 *
 * <h3>General Information</h3>
 * <p>
 * This tool can download and package data sources for a user-specified species.
 * The data sources downloaded by this tool correspond to the latest /current versions of the data sources supported as defined in Ensembl database.
 * </p>
 *
 * <p>
 * To download, package and extract the data sources, you can invoke {@link FuncotatorDataSourceBundler} in the following way:
 *      <ul>
 *          ./gatk FuncotatorDataSourceBundler \
 *          -N organismName \
 *          -S subgroup \
 *          -O outputFile \
 *          --overwrite-output-file
 *          --extract-data-source
 *      </ul>
 * </p>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>By default {@link FuncotatorDataSourceBundler} will not overwrite data sources if they already exist locally. </li>
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Download and package data sources for a given organism to be used for Funcotator.",
        oneLineSummary = "Data source bundler for Funcotator.",
        programGroup = VariantEvaluationProgramGroup.class //Need to make a new program group
)
@DocumentedFeature
public class FuncotatorDataSourceBundler extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(FuncotatorDataSourceBundler.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String OVERWRITE_ARG_LONG_NAME = "overwrite-output-file";
    public static final String EXTRACT_AFTER_DOWNLOAD  = "extract-after-download";

    //==================================================================================================================
    // Private Static Members:

    // Set to always get the latest version of the data sources:
    private static final String BASE_URL = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION;

    //will maybe add in variables for the urls for each of the different organisms

    //==================================================================================================================
    // Private Members:

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName  = OVERWRITE_ARG_LONG_NAME,
            doc = "Overwrite output file if it exists already.",
            optional = true)
    private boolean overwriteOutputFile = false;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output location for the data sources.",
            optional = true)
    protected File outputFile;

    @Argument(
            shortName = EXTRACT_AFTER_DOWNLOAD,
            fullName  = EXTRACT_AFTER_DOWNLOAD,
            doc = "Extract the data sources to a sibling folder after they have been downloaded.",
            optional = true)
    protected boolean extractDataSourcesAfterDownload = false;

    @Argument(
            shortName = StandardArgumentDefinitions.ORGANISM_TYPE_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.ORGANISM_TYPE_LONG_NAME,
            doc = "Organism we want to download data sources for.")
    protected String organismName;

    @Argument(
            shortName = StandardArgumentDefinitions.SUBGROUP_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.SUBGROUP_LONG_NAME,
            doc = "Subgroup of organism type which we want the data sources for.")
    protected String subgroupName;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    protected void onStartup() {

        // Make sure the user specified an organism type and a subgroup type for the data source:
        if ((organismName == null) || (subgroupName == null)) {
            throw new UserException("Must select an organism and subgroup for data source.");
        }

        // Make sure the testing inputs are correct:
        if ( ((organismName == null) && (subgroupName != null)) || ((subgroupName == null) && (organismName != null)) ) {
            throw new UserException("Must specify both an organism type and a subgroup type.");
        }

        if ( overwriteOutputFile ) {
            logger.info("Overwrite ENABLED. Will overwrite existing data sources download.");
        }
    }

    @Override
    protected Object doWork() {

        final String dataSourcePath = null;
        return dataSourcePath;
    }

}
