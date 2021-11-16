package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import picard.sam.CreateSequenceDictionary;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Calendar;


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
 *          -organismName \
 *          -species-name speciesName \
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
@BetaFeature
public class FuncotatorDataSourceBundler extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(FuncotatorDataSourceBundler.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String BACTERIA_ARG_LONG_NAME   = "bacteria";
    public static final String FUNGI_ARG_LONG_NAME      = "fungi";
    public static final String METAZOA_ARG_LONG_NAME    = "metazoa";
    public static final String PLANTS_ARG_LONG_NAME     = "plants";
    public static final String PROTISTS_ARG_LONG_NAME   = "protists";
    public static final String SPECIES_ARG_LONG_NAME    = "species-name";
    public static final String OVERWRITE_ARG_LONG_NAME            = "overwrite-output-file";
    public static final String OUTPUT_DATASOURCES_FOLDER_ARG_NAME = "output-datasources-folder";

    //==================================================================================================================
    // Private Static Members:

    // Set to always get the latest version of the data sources:
    private static final String  BASE_URL = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION;

    // Will maybe make these private again
    @VisibleForTesting
    public static final String BACTERIA_BASE_URL      = BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION;
    public static final String FUNGI_BASE_URL         = BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;
    public static final String METAZOA_BASE_URL       = BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;
    public static final String PLANTS_BASE_URL        = BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;
    public static final String PROTISTS_BASE_URL      = BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;

    public static final String BACTERIA_BASE_FASTA    = BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION;
    public static final String FUNGI_BASE_FASTA       = BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION;
    public static final String METAZOA_BASE_FASTA     = BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION;
    public static final String PLANTS_BASE_FASTA      = BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION;
    public static final String PROTISTS_BASE_FASTA    = BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION;

    //==================================================================================================================
    // Private Static Members:
    protected static final int BUFFER_SIZE_BYTES    = 1024 * 1024;


    //==================================================================================================================
    // Private Members:

    // Arguments:
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
            shortName = SPECIES_ARG_LONG_NAME,
            fullName  = SPECIES_ARG_LONG_NAME,
            doc = "Download data sources for this species of the organism.")
    protected String speciesName;

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName  = OVERWRITE_ARG_LONG_NAME,
            doc = "Overwrite output file if it exists already.",
            optional = true)
    private boolean overwriteOutputFile = false;

    @Argument(fullName = OUTPUT_DATASOURCES_FOLDER_ARG_NAME,
            shortName  = OUTPUT_DATASOURCES_FOLDER_ARG_NAME,
            doc = "Location in which to put the output datasources.",
            optional = true)
    private String outputDatasourcesFolder = null;

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
            throw new UserException("Must specify a species to bundle data sources for.");
        }

        if (overwriteOutputFile) {
            logger.info("Overwrite ENABLED. Will overwrite existing data sources download.");
        }

    }

    @Override
    protected Object doWork() {

        final String dataSourceOrganism;
        final String baseURL;
        final String baseFastaURL;

        // Get the correct data source:
        if ( getBacteriaDataSources ) {
            dataSourceOrganism = "bacteria";
            baseURL = BACTERIA_BASE_URL;
            baseFastaURL = BACTERIA_BASE_FASTA;

        } else if ( getFungiDataSources ) {
            dataSourceOrganism = "fungi";
            baseURL = FUNGI_BASE_URL;
            baseFastaURL = FUNGI_BASE_FASTA;

        } else if ( getMetazoaDataSources ) {
            dataSourceOrganism = "metazoa";
            baseURL = METAZOA_BASE_URL;
            baseFastaURL = METAZOA_BASE_FASTA;

        } else if ( getPlantsDataSources ) {
            dataSourceOrganism = "plants";
            baseURL = PLANTS_BASE_URL;
            baseFastaURL = PLANTS_BASE_FASTA;

        } else {
            dataSourceOrganism = "protists";
            baseURL = PROTISTS_BASE_URL;
            baseFastaURL = PROTISTS_BASE_FASTA;
        }

        downloadAndValidateDataSources(dataSourceOrganism, speciesName, baseURL, baseFastaURL);

        // Token return value:
        return true;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private void downloadAndValidateDataSources(final String dsOrganism, final String dsSpecies, final String baseURL, final String baseFastaURL) {

        logger.info(dsOrganism + ":" + dsSpecies + " data sources selected.");

        // Get or create a named location to put our data sources:
        final String outputFolder;
        if (outputDatasourcesFolder != null) {
            outputFolder = outputDatasourcesFolder;
        }
        else {
            outputFolder = speciesName + "_dataSources.v0.0." + FuncotatorDataSourceBundlerUtils.getCurrentDateString();
        }

        // Make folders to put data sources in:
        logger.info("Creating output folder for datasources: " + outputFolder);
        makeDataSourcesFolderStructure(outputFolder, speciesName);

        // Make the bundler object:
        final FuncotatorDataSourceBundlerHttpClient bundler = new FuncotatorDataSourceBundlerHttpClient(
                IOUtils.getPath(outputFolder), dsOrganism, speciesName, baseURL, baseFastaURL
        );

        // ===================================
        // Download the data sources:
        logger.info("Downloading data files...");
        bundler.downloadDataSources();

        // ===================================
        // Process downloaded files:

        // Extract our data:
        logger.info("Extracting gzipped data...");
        bundler.extractGzippedFiles(overwriteOutputFile);

        // Build the fasta index file:
        logger.info("Indexing FASTA file...");
        try {
            FastaSequenceIndexCreator.create(bundler.getFastaUnzipPath(), true);
        } catch (IOException e) {
            throw new UserException("Error. Unable to create index for fasta file.", e);
        }

        // Build the fasta .dict file:
        logger.info("Creating sequence dictionary...");
        final CreateSequenceDictionary csd_tool = new CreateSequenceDictionary();
        csd_tool.instanceMain(new String[] {"-R", bundler.getFastaUnzipPath().toString()});

        // Index the gtf file:
        logger.info("Indexing GTF file...");
        bundler.sortAndIndexGtfFile();

        // Build the config file for the new data source directory:
        logger.info("Building GTF data source config file...");
        bundler.buildConfigFile();

        // Create a manifest file:
        logger.info("Building datasources manifest file...");
        bundler.buildManifestFile();

        // Create the template config file:
        logger.info("Creating template config file...");
        bundler.buildTemplateConfigFile();

        // Create README file:
        logger.info("Creating high-level README file...");
        bundler.buildReadMeFile();

        // ===================================
        // CLEANUP:
        logger.info("Cleaning up intermediate files...");

        // Delete gtf and fasta zip files
        IOUtils.deleteOnExit(bundler.getOutputDestination());
        IOUtils.deleteOnExit(bundler.getFastaOutputDestination());

        // Delete gtf file:
        IOUtils.deleteOnExit(bundler.getDSUnzipPath());
    }

    /**
     * Create a folder with the given path or raise an exception.
     * @param baseFolder The name of the base folder for the new data sources.
     */
    private static void makeFolderOrFail(final Path baseFolder) {
        final boolean success = baseFolder.toAbsolutePath().toFile().mkdir();
        if (!success) {
            throw new UserException("Unable to make file.");
        }
    }

    /**
     * Method which builds the folder structure for the funcotator data sources.
     * @param baseFolder The name of the base folder for the new data sources.
     * @param speciesName The name of the species for which to create data sources.
     */
    public static void makeDataSourcesFolderStructure(final String baseFolder, final String speciesName) {
        // Base folder:
        makeFolderOrFail(IOUtils.getPath(baseFolder));

        // Subfolder for GTF:
        makeFolderOrFail(IOUtils.getPath(baseFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION));

        // Subfolder for species:
        makeFolderOrFail(IOUtils.getPath(baseFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName));
    }
}



