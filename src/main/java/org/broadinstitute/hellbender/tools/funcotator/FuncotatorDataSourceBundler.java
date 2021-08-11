package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
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
    public static final String OVERWRITE_ARG_LONG_NAME  = "overwrite-output-file";
    public static final String EXTRACT_AFTER_DOWNLOAD   = "extract-after-download";

    //==================================================================================================================
    // Private Static Members:

    // Set to always get the latest version of the data sources:
    private static final String BASE_URL = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION;

    // Will maybe make these private again
    @VisibleForTesting
    public static final String BACTERIA_BASE_URL    = BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION;
    public static final String FUNGI_BASE_URL       = BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;
    public static final String METAZOA_BASE_URL     = BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;
    public static final String PLANTS_BASE_URL      = BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;
    public static final String PROTISTS_BASE_URL    = BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;

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

    // Copy buffer:
    protected final byte copyBuffer[] = new byte[BUFFER_SIZE_BYTES];

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
            fullName = SPECIES_ARG_LONG_NAME,
            doc = "Download data sources for this species of the organism.")
    protected String speciesName;

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName = OVERWRITE_ARG_LONG_NAME,
            doc = "Overwrite output file if it exists already.",
            optional = true)
    private boolean overwriteOutputFile = false;

    @Argument(
            shortName = EXTRACT_AFTER_DOWNLOAD,
            fullName = EXTRACT_AFTER_DOWNLOAD,
            doc = "Extract the data sources to a sibling folder after they have been downloaded.",
            optional = true)
    protected boolean extractDataSourcesAfterDownload = false;

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
        if (getBacteriaDataSources) {
            dataSourceOrganism = "bacteria";
            baseURL = BACTERIA_BASE_URL;
            baseFastaURL = BACTERIA_BASE_FASTA;
        } else if (getFungiDataSources) {
            dataSourceOrganism = "fungi";
            baseURL = FUNGI_BASE_URL;
            baseFastaURL = FUNGI_BASE_FASTA;
        } else if (getMetazoaDataSources) {
            dataSourceOrganism = "metazoa";
            baseURL = METAZOA_BASE_URL;
            baseFastaURL = METAZOA_BASE_FASTA;
        } else if (getPlantsDataSources) {
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

        // Make folders to put data sources in:
        makeFolders(IOUtils.getPath(speciesName));

        // Make the bundler object:
        final FuncotatorDataSourceBundlerHttpClient bundler = FuncotatorDataSourceBundlerHttpClient.create(dsOrganism, speciesName, baseURL, baseFastaURL);


        // Download the gtf file:
        FuncotatorDataSourceBundlerHttpClient.downloadDataSources(bundler.getDSUrl(), bundler.getOutputDestination());

        // Download the fasta file:
        FuncotatorDataSourceBundlerHttpClient.downloadDataSources(bundler.getFastaURL(), bundler.getFastaOutputDestination());


        // Extract data sources if requested:
        if ( extractDataSourcesAfterDownload ) {
            FuncotatorDataSourceBundlerUtils.extractGtfGz(bundler.getOutputDestination().toString(), bundler.getDSUnzipPath().toString(), overwriteOutputFile);
            FuncotatorDataSourceBundlerUtils.extractGtfGz(bundler.getFastaOutputDestination().toString(), bundler.getFastaUnzipPath().toString(), overwriteOutputFile);
        }
        else {
            logger.info("IMPORTANT: You must unzip the downloaded data sources prior to using them with Funcotator.");
        }


        // Index the fasta file and build the dict file:
        FuncotatorDataSourceBundlerHttpClient.buildFastaIndexFile(bundler.getFastaUnzipPath());

        // Copy the config file into new directory:
        FuncotatorDataSourceBundlerHttpClient.buildConfigFile(bundler);

        // Create a manifest file:
        FuncotatorDataSourceBundlerHttpClient.buildManifestFile(bundler);

        // Create the template config file:
        FuncotatorDataSourceBundlerHttpClient.buildTemplateConfigFile(bundler);

        // Create ReadMe file:
        FuncotatorDataSourceBundlerHttpClient.buildReadMeFile(bundler);

        // Download the gtf ReadMe file for specific data source file:
        FuncotatorDataSourceBundlerHttpClient.downloadDataSources(bundler.getGtfReadMeURL(), bundler.getGtfReadMePath());

        // Download the fasta ReadMe file for specific data source file:
        FuncotatorDataSourceBundlerHttpClient.downloadDataSources(bundler.getFastaReadMeURL(), bundler.getFastaReadMePath());

//        // Index the gtf file:
//        FuncotatorDataSourceBundlerHttpClient.buildIndexFile(bundler.getDSUnzipPath(), bundler.getIndexPath());

    }

    public static void makeFolders(Path speciesName) {
        String path = speciesName.toAbsolutePath().toString();
        File newFolder = new File(path);
        boolean bool = newFolder.mkdir();
        if (!bool) {
            throw new UserException("Unable to make file.");
        }
        makeFolder2(path, speciesName);
    }

    public static void makeFolder2(String pathName, Path speciesName) {
        String path = pathName + "/" + DataSourceUtils.ENSEMBL_EXTENSION;
        File newFolder = new File(path);
        boolean bool = newFolder.mkdir();
        if (!bool) {
            throw new UserException("Unable to make file.");
        }
        makeFolder3(path, speciesName);
    }

    public static void makeFolder3(String pathName, Path speciesName) {
        String path = pathName + "/" + speciesName.toString();
        File newFolder = new File(path);
        boolean bool = newFolder.mkdir();
        if (!bool) {
            throw new UserException("Unable to make file.");
        }
    }


    //==================================================================================================================
    // Helper Data Types:
}



