package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import picard.sam.CreateSequenceDictionary;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;


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
 * The currently supported kingdoms are:
 * </p>
 * <ul>
 *     <li>bacteria</li>
 *     <li>protists</li>
 *     <li>fungi</li>
 *     <li>plants</li>
 *     <li>metazoa</li>
 * </ul>
 *
 * <p>
 *     A list of supported species for each kingdom can be found at the following web addresses:
 * </p>
 * <table>
 *     <th>
 *         <td>Kingdom</td><td>Web Address(es)</td>
 *     </th>
 *     <tr>
 *         <td>bacteria</td>
 *         <td><a href="http://ftp.ensemblgenomes.org/pub/bacteria/current/uniprot_report_EnsemblBacteria.txt">http://ftp.ensemblgenomes.org/pub/bacteria/current/uniprot_report_EnsemblBacteria.txt</a></td>
 *     </tr>
 *     <tr>
 *         <td>protists</td>
 *         <td><a href="http://ftp.ensemblgenomes.org/pub/protists/current/uniprot_report_EnsemblProtists.txt">http://ftp.ensemblgenomes.org/pub/protists/current/uniprot_report_EnsemblProtists.txt</a></td>
 *     </tr>
 *     <tr>
 *         <td>fungi</td>
 *         <td><a href="http://ftp.ensemblgenomes.org/pub/fungi/current/uniprot_report_EnsemblFungi.txt">http://ftp.ensemblgenomes.org/pub/fungi/current/uniprot_report_EnsemblFungi.txt</a></td>
 *     </tr>
 *     <tr>
 *         <td>plants</td>
 *         <td><a href="http://ftp.ensemblgenomes.org/pub/plants/current/uniprot_report_EnsemblPlants.txt">http://ftp.ensemblgenomes.org/pub/plants/current/uniprot_report_EnsemblPlants.txt</a></td>
 *     </tr>
 *     <tr>
 *         <td>metazoa</td>
 *         <td><a href="http://ftp.ensemblgenomes.org/pub/metazoa/current/uniprot_report_EnsemblMetazoa.txt">http://ftp.ensemblgenomes.org/pub/metazoa/current/uniprot_report_EnsemblMetazoa.txt</a></td>
 *     </tr>
 * </table>
 *
 * <p>
 * To download, package and extract the data sources, you can invoke {@link FuncotatorDataSourceBundler} as follows:
 *      <ul>
 *          ./gatk FuncotatorDataSourceBundler \
 *          --organism-kingdom kingdom \
 *          --species-name speciesName \
 *          -O outputFile \
 *          --overwrite-output-file
 *          --extract-data-source
 *      </ul>
 * Note: Both a kingdom and a species are required.
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

    public static final String ORGANISM_KINGDOM_ARG_LONG_NAME = "organism-kingdom";
    public static final String SPECIES_ARG_LONG_NAME          = "species-name";
    public static final String OVERWRITE_ARG_LONG_NAME     = "overwrite-output-file";

    // Set to always get the latest version of the data sources:
    private static final String BASE_URL = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION;

    public enum OrganismKingdom {
        BACTERIA(BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION,
                BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION),
        FUNGI(BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION,
                BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION) ,
        METAZOA(BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION,
                BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION),
        PLANTS(BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION,
                BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION),
        PROTISTS(BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION,
                BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION + DataSourceUtils.FASTA_EXTENSION);

        OrganismKingdom(final String baseUrl, final String baseFasta) {
            this.baseUrl = baseUrl;
            this.baseFastaUrl = baseFasta;
            uniprotFileName = "uniprot_report_Ensembl" + this.name().substring(0,1).toUpperCase() + this.name().substring(1).toLowerCase() + ".txt";
        }

        private final String baseUrl;
        private final String baseFastaUrl;
        private final String uniprotFileName;

        public String getBaseUrl() { return baseUrl; }
        public String getBaseFastaUrl() { return baseFastaUrl; }
        public String getUniprotFileName() { return uniprotFileName; }

        @Override
        public String toString() {
            return this.name().toLowerCase();
        }
    }

    //==================================================================================================================
    // Private Members:

    // Arguments:
    @Argument(fullName = ORGANISM_KINGDOM_ARG_LONG_NAME,
            shortName = ORGANISM_KINGDOM_ARG_LONG_NAME,
            doc = "Kingdom of the organism for which to download datasources.")
    private OrganismKingdom kingdom;

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

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output location for the data sources.  Cannot be a bucket.",
            optional = true)
    private String outputDatasourcesFolder = null;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    protected void onStartup() {
        if (overwriteOutputFile) {
            logger.info("Overwrite ENABLED. Will overwrite existing data sources download.");
        }
    }

    @Override
    protected Object doWork() {
        downloadAndValidateDataSources();

        // Token return value:
        return true;
    }

    //==================================================================================================================
    // Instance Methods:

    private void downloadAndValidateDataSources() {

        logger.info(kingdom.toString() + ":" + speciesName + " data sources selected.");

        // Get or create a named location to put our data sources:
        final String outputFolder;
        if (outputDatasourcesFolder != null) {
            outputFolder = outputDatasourcesFolder;
        }
        else {
            outputFolder = this.speciesName + "_dataSources.v0.0." + FuncotatorDataSourceBundlerUtils.getCurrentDateString();
        }

        // Check if the output folder exists:
        if (Files.exists(IOUtils.getPath(outputFolder)) && !overwriteOutputFile) {
            logger.error("Output folder already exists: " + outputFolder);
            throw new UserException("Output folder for datasources already exists (" + outputFolder + ").  Cannot continue.");
        }

        // Make folders to put data sources in:
        logger.info("Creating output folder for datasources: " + outputFolder);
        makeDataSourcesFolderStructure(outputFolder, this.speciesName, overwriteOutputFile);

        // Make the bundler object:
        final FuncotatorDataSourceBundlerHttpClient bundler = new FuncotatorDataSourceBundlerHttpClient(
                IOUtils.getPath(outputFolder), kingdom, this.speciesName
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
        try {
            Files.delete(bundler.getOutputDestination());
            Files.delete(bundler.getFastaOutputDestination());

            // Delete gtf file:
            Files.delete(bundler.getDSUnzipPath());
        } catch (final IOException ex) {
            throw new GATKException("ERROR: Could not clean up intermediate files!", ex);
        }
    }

    /**
     * Create a folder with the given path or raise an exception.
     * @param baseFolder The name of the base folder for the new data sources.
     */
    private static void makeFolderOrFail(final Path baseFolder, final Boolean existOK) {

        if (Files.exists(baseFolder) && existOK) {
            return;
        }

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
    public static void makeDataSourcesFolderStructure(final String baseFolder, final String speciesName, final Boolean overwriteOutputFile) {
        // Base folder:
        makeFolderOrFail(IOUtils.getPath(baseFolder), overwriteOutputFile);

        // Subfolder for GTF:
        makeFolderOrFail(IOUtils.getPath(baseFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION), overwriteOutputFile);

        // Subfolder for species:
        makeFolderOrFail(IOUtils.getPath(baseFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName), overwriteOutputFile);
    }
}



