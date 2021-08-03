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

    public static final String BACTERIA_ARG_LONG_NAME           = "bacteria";
    public static final String FUNGI_ARG_LONG_NAME              = "fungi";
    public static final String METAZOA_ARG_LONG_NAME            = "metazoa";
    public static final String PLANTS_ARG_LONG_NAME             = "plants";
    public static final String PROTISTS_ARG_LONG_NAME           = "protists";
    public static final String OVERWRITE_ARG_LONG_NAME          = "overwrite-output-file";
    public static final String EXTRACT_AFTER_DOWNLOAD           = "extract-after-download";
    public static final String VALIDATE_INTEGRITY_ARG_LONG_NAME = "validate-integrity";

    //==================================================================================================================
    // Private Static Members:

    // Set to always get the latest version of the data sources:
    private static final String BASE_URL = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION;

    private static final String BACTERIA_BASE_URL = BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION;

    private static final Path BACTERIA_PATH = IOUtils.getPath(BACTERIA_BASE_URL + DataSourceUtils.GTF_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION);

    private static final String FUNGI_BASE_URL = BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION;

    private static final Path FUNGI_PATH = IOUtils.getPath(FUNGI_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    private static final String METAZOA_BASE_URL = BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION;

    private static final Path METAZOA_PATH = IOUtils.getPath(METAZOA_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    private static final String PLANTS_BASE_URL = BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION;

    private static final Path PLANTS_PATH = IOUtils.getPath(PLANTS_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    private static final String PROTISTS_BASE_URL = BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION;

    private static final Path PROTISTS_PATH = IOUtils.getPath(PROTISTS_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    //will maybe add in variables for the urls for each of the different organisms

    //==================================================================================================================
    // Private Members:

    @Argument(fullName = BACTERIA_ARG_LONG_NAME,
            shortName  = BACTERIA_ARG_LONG_NAME,
            doc = "Download data sources for bacteria.",
            optional = true)
    private boolean getBacteriaDataSources = false;

    @Argument(fullName = FUNGI_ARG_LONG_NAME,
            shortName  = FUNGI_ARG_LONG_NAME,
            doc = "Download data sources for fungi.",
            optional = true)
    private boolean getFungiDataSources = false;

    @Argument(fullName = METAZOA_ARG_LONG_NAME,
            shortName  = METAZOA_ARG_LONG_NAME,
            doc = "Download data sources for metazoa.",
            optional = true)
    private boolean getMetazoaDataSources = false;

    @Argument(fullName = PLANTS_ARG_LONG_NAME,
            shortName  = PLANTS_ARG_LONG_NAME,
            doc = "Download data sources for plants.",
            optional = true)
    private boolean getPlantsDataSources = false;

    @Argument(fullName = PROTISTS_ARG_LONG_NAME,
            shortName  = PROTISTS_ARG_LONG_NAME,
            doc = "Download data sources for protists.",
            optional = true)
    private boolean getProtistsDataSources = false;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output location for the data sources.",
            optional = true)
    protected File outputFile;

    @Argument(
            shortName = StandardArgumentDefinitions.SPECIES_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.SPECIES_LONG_NAME,
            doc = "Download data sources for this species of the organism.")
    protected String speciesName;

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName  = OVERWRITE_ARG_LONG_NAME,
            doc = "Overwrite output file if it exists already.",
            optional = true)
    private boolean overwriteOutputFile = false;

    @Argument(
            shortName = EXTRACT_AFTER_DOWNLOAD,
            fullName  = EXTRACT_AFTER_DOWNLOAD,
            doc = "Extract the data sources to a sibling folder after they have been downloaded.",
            optional = true)
    protected boolean extractDataSourcesAfterDownload = false;

    @Argument(fullName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            shortName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            doc = "Validate the integrity of the data sources after downloading them using sha256.",
            optional = true)
    private boolean doValidateIntegrity = false;

//    @Argument(
//            shortName = StandardArgumentDefinitions.ORGANISM_TYPE_SHORT_NAME,
//            fullName  = StandardArgumentDefinitions.ORGANISM_TYPE_LONG_NAME,
//            doc = "Organism we want to download data sources for.")
//    protected String organismName;
//
//    @Argument(
//            shortName = StandardArgumentDefinitions.SUBGROUP_SHORT_NAME,
//            fullName  = StandardArgumentDefinitions.SUBGROUP_LONG_NAME,
//            doc = "Subgroup of organism type which we want the data sources for.")
//    protected String subgroupName;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    protected void onStartup() {

        // Make sure the user specified an organism data source to bundle
        if ((!getBacteriaDataSources) && (!getFungiDataSources) && (!getMetazoaDataSources) && (!getPlantsDataSources) && (!getProtistsDataSources)) {
            throw new UserException("Must select either bacteria, fungi, metazoa, plants, or protists data sources.");
        }

        // Make sure the user specified a species data source to bundle
        if (speciesName == null) {
            throw new UserException("Must specify a species to bundle data sources for."):
        }

//        // Make sure the user specified an organism type and a subgroup type for the data source:
//        if ((organismName == null) || (subgroupName == null)) {
//            throw new UserException("Must select an organism and subgroup for data source.");
//        }
//
//        // Make sure the testing inputs are correct:
//        if ( ((organismName == null) && (subgroupName != null)) || ((subgroupName == null) && (organismName != null)) ) {
//            throw new UserException("Must specify both an organism type and a subgroup type.");
//        }

        if ( overwriteOutputFile ) {
            logger.info("Overwrite ENABLED. Will overwrite existing data sources download.");
        }
    }

    @Override
    protected Object doWork() {

        final String dataSourceOrganism;
        final Path dataSourcePath;

        // Get the correct data source:
        if ( getBacteriaDataSources ) {
            dataSourceOrganism = "Bacteria";
            dataSourcePath = BACTERIA_PATH;
        }
        else if ( getFungiDataSources ) {
            dataSourceOrganism = "Fungi";
            dataSourcePath = FUNGI_PATH;
        }
        else if ( getMetazoaDataSources ) {
            dataSourceOrganism = "Metazoa";
            dataSourcePath = METAZOA_PATH;
        }
        else if ( getPlantsDataSources ) {
            dataSourceOrganism = "Plants";
            dataSourcePath = PLANTS_PATH;
        }
        else {
            dataSourceOrganism = "Protists";
            dataSourcePath = PROTISTS_PATH;
        }

        downloadAndValidateDataSources(dataSourceOrganism, speciesName, dataSourcePath);

        // Token return value:
        return true;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private NioFileCopierWithProgressMeter createNioDownloader(final Path dsPath) {
        // Get the data sources file:
        final Path outputDestination = getOutputLocation(dsPath);

        // Set up and initiate our download:
        return NioFileCopierWithProgressMeter.create(dsPath, outputDestination, overwriteOutputFile);

    }

    private void downloadAndValidateDataSources(final String dataSourceType, final String speciesName, final Path dsPath) {
        logger.info(dataSourceType + " data sources selected.");
        // Confirm file integrity if requested:
//        if ( doValidateIntegrity ) {
//            validateIntegrity(results);
//        }
        // Create our downloader:
        final NioFileCopierWithProgressMeter xerox = createNioDownloader(dsPath);

        // Get the datasources file:
        final NioFileCopierWithProgressMeterResults results = downloadDataSources(xerox);

//

        // Extract data sources if requested:
        if ( extractDataSourcesAfterDownload ) {
            IOUtils.extractTarGz(results.getDestination(), results.getDestination().getParent(), overwriteOutputFile);
        }
        else {
            logger.info("IMPORTANT: You must unzip the downloaded data sources prior to using them with Funcotator.");
        }
    }

    private NioFileCopierWithProgressMeterResults downloadDataSources(final NioFileCopierWithProgressMeter xerox) {

        // Set up the validity check while the download occurs if requested:
//        if ( doValidateIntegrity ) {
//            // Read the sha256sum into memory:
//            final String expectedSha256Sum = readSha256SumFromPath(checksumPath);
//
//            // Setup the copier to calculate the checksum:
//            xerox.setChecksumAlgorithmAndExpectedChecksum(MessageDigestAlgorithms.SHA_256, expectedSha256Sum);
//        }

        // Initiate the copy and return our results:
        return xerox.initiateCopy();

    }

    private Path getOutputLocation(final Path dataSourcesPath) {
        if ( outputFile == null ) {
            return IOUtils.getPath(dataSourcesPath.getFileName().toString());
        }
        else {
            return outputFile.toPath();
        }
    }

//        private void validateIntegrity(final NioFileCopierWithProgressMeterResults results) {
//
//            // verify the hashes are the same:
//            if ( !results.isDestFileValid() ) {
//                throw new UserException("ERROR: downloaded data sources are corrupt!  Unexpected checksum: " + results.getChecksum() + " != " + results.getExpectedChecksum());
//            }
//            else {
//                logger.info("Integrity check on downloaded data sources succeeded.");
//            }
//        }

//        private String readSha256SumFromPath(final Path sha256Path) {
//            final String expectedSha256Sum;
//            try {
//                logger.info("Collecting expected checksum...");
//                expectedSha256Sum = Files.lines(sha256Path).findFirst().orElse(null);
//
//                if ( expectedSha256Sum == null ) {
//                    throw new UserException("Unable to retrieve expected checksum from: " + sha256Path.toUri());
//                }
//
//                logger.info("Collection complete!");
//            }
//            catch ( final IOException ex ) {
//                throw new UserException("Could not read in sha256sum from file: " + sha256Path.toUri().toString(), ex);
//            }
//
//            // Clean up and return the checksum:
//            return cleanExpectedSha256SumString(expectedSha256Sum);
        }


//        private String cleanExpectedSha256SumString(final String expectedSha256SumString) {
//
//            String cleanString = expectedSha256SumString.trim().toLowerCase();
//
//            // The format of the file can contain the checksum, followed by the file name.
//            // If this is the case, we need to truncate the string:
//            if ( cleanString.contains(" ") ) {
//                cleanString = cleanString.substring(0, cleanString.indexOf(" "));
//            }
//            if ( cleanString.contains("\t")) {
//                cleanString = cleanString.substring(0, cleanString.indexOf("\t"));
//            }
//
//            return cleanString;
//        }

        //==================================================================================================================
        // Helper Data Types:



