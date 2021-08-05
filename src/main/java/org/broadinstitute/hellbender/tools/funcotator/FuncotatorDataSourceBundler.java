package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import org.apache.arrow.flatbuf.Int;
import org.apache.http.HttpClientConnection;
//import org.apache.commons.httpclient.*;
//import org.apache.commons.httpclient.methods.*;
import sun.net.www.http.HttpClient;
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
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeter;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeterResults;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.checkerframework.common.reflection.qual.GetMethod;
import com.google.cloud.storage.HttpMethod;
//import org.apache.commons.httpclient.methods;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.HttpEntity;
import org.apache.http.NameValuePair;
import org.apache.http.client.entity.UrlEncodedFormEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.message.BasicNameValuePair;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpUriRequest;

import java.io.FileOutputStream;
import java.io.IOException;
import org.apache.http.client.ClientProtocolException;

import java.io.InputStream;
import java.io.OutputStream;

import java.net.URL;
import java.net.HttpURLConnection;

import java.io.File;
import java.nio.file.Path;

//import org.jsoup.Connection;
//import org.jsoup.Jsoup;
//import org.jsoup.nodes.Document;

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

    public static final String BACTERIA_ARG_LONG_NAME = "bacteria";
    public static final String FUNGI_ARG_LONG_NAME = "fungi";
    public static final String METAZOA_ARG_LONG_NAME = "metazoa";
    public static final String PLANTS_ARG_LONG_NAME = "plants";
    public static final String PROTISTS_ARG_LONG_NAME = "protists";
    public static final String SPECIES_ARG_LONG_NAME = "species-name";
    public static final String OVERWRITE_ARG_LONG_NAME = "overwrite-output-file";
    public static final String EXTRACT_AFTER_DOWNLOAD = "extract-after-download";

    //==================================================================================================================
    // Private Static Members:

    // Set to always get the latest version of the data sources:
    private static final String BASE_URL = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION;

    private static final String BACTERIA_BASE_URL = BASE_URL + DataSourceUtils.BACTERIA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION;

//    private static final Path BACTERIA_PATH = IOUtils.getPath(BACTERIA_BASE_URL + DataSourceUtils.GTF_EXTENSION + DataSourceUtils.BACTERIA_COLLECTION_EXTENSION);

    private static final String FUNGI_BASE_URL = BASE_URL + DataSourceUtils.FUNGI_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;

//    private static final Path FUNGI_PATH = IOUtils.getPath(FUNGI_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    private static final String METAZOA_BASE_URL = BASE_URL + DataSourceUtils.METAZOA_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;

//    private static final Path METAZOA_PATH = IOUtils.getPath(METAZOA_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    private static final String PLANTS_BASE_URL = BASE_URL + DataSourceUtils.PLANTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;

//    private static final Path PLANTS_PATH = IOUtils.getPath(PLANTS_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    private static final String PROTISTS_BASE_URL = BASE_URL + DataSourceUtils.PROTISTS_DS_EXTENSION + DataSourceUtils.GTF_EXTENSION;

//    private static final Path PROTISTS_PATH = IOUtils.getPath(PROTISTS_BASE_URL + DataSourceUtils.GTF_EXTENSION);

    protected static final int BUFFER_SIZE_BYTES = 1024 * 1024;

    //will maybe add in variables for the urls for each of the different organisms

    //==================================================================================================================
    // Private Members:

    // Copy buffer:
    protected final byte copyBuffer[] = new byte[BUFFER_SIZE_BYTES];

    // Arguments:
    @Argument(fullName = BACTERIA_ARG_LONG_NAME,
            shortName = BACTERIA_ARG_LONG_NAME,
            doc = "Download data sources for bacteria.",
            optional = true)
    private boolean getBacteriaDataSources = false;

    @Argument(fullName = FUNGI_ARG_LONG_NAME,
            shortName = FUNGI_ARG_LONG_NAME,
            doc = "Download data sources for fungi.",
            optional = true)
    private boolean getFungiDataSources = false;

    @Argument(fullName = METAZOA_ARG_LONG_NAME,
            shortName = METAZOA_ARG_LONG_NAME,
            doc = "Download data sources for metazoa.",
            optional = true)
    private boolean getMetazoaDataSources = false;

    @Argument(fullName = PLANTS_ARG_LONG_NAME,
            shortName = PLANTS_ARG_LONG_NAME,
            doc = "Download data sources for plants.",
            optional = true)
    private boolean getPlantsDataSources = false;

    @Argument(fullName = PROTISTS_ARG_LONG_NAME,
            shortName = PROTISTS_ARG_LONG_NAME,
            doc = "Download data sources for protists.",
            optional = true)
    private boolean getProtistsDataSources = false;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
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
        final String dataSourceSpecies = speciesName;
        final String dataSourceURL;
        final Path dataSourcePath;

        // Get the correct data source:
        if (getBacteriaDataSources) {
            dataSourceOrganism = "bcteria";
            dataSourceURL = BACTERIA_BASE_URL + speciesName + "/";
            dataSourcePath = IOUtils.getPath(BACTERIA_BASE_URL + speciesName "/" );
        } else if (getFungiDataSources) {
            dataSourceOrganism = "fungi";
            dataSourceURL = FUNGI_BASE_URL + speciesName + "/";
            dataSourcePath = IOUtils.getPath(FUNGI_BASE_URL + speciesName + "/");
        } else if (getMetazoaDataSources) {
            dataSourceOrganism = "metazoa";
            dataSourceURL = METAZOA_BASE_URL + speciesName + "/";
            dataSourcePath = IOUtils.getPath(METAZOA_BASE_URL + speciesName + "/");
        } else if (getPlantsDataSources) {
            dataSourceOrganism = "plants";
            dataSourceURL = PLANTS_BASE_URL + speciesName + "/";
            dataSourcePath = IOUtils.getPath(PLANTS_BASE_URL + speciesName + "/");
        } else {
            dataSourceOrganism = "protists";
            dataSourceURL = PROTISTS_BASE_URL + speciesName + "/";
            dataSourcePath = IOUtils.getPath(PROTISTS_BASE_URL + speciesName + "/");
        }

        downloadAndValidateDataSources(dataSourceOrganism, dataSourceSpecies, dataSourceURL, dataSourcePath);

        // Token return value:
        return true;
    }

    //==================================================================================================================
    // Static Methods:

    @VisibleForTesting
    static Path getPath(String organismBaseURL, String speciesName) {
        Path testPath = IOUtils.getPath(organismBaseURL + speciesName + "/");
        return testPath;
    }

    //==================================================================================================================
    // Instance Methods:

    private void downloadAndValidateDataSources(final String dsOrganism, final String dsSpecies, final String dsURL, final Path dsPath){
        logger.info(dsOrganism + ":" + dsSpecies + " data sources selected.");

        // Creating CloseableHttpClient object to access the webpage and retrieve the file:
        CloseableHttpClient client = HttpClientBuilder.create().build();

        // Creating an HttpGet object to send the request to the server:
        HttpGet request = new HttpGet(dsURL + "/" + dsSpecies); //instead of dsSpecies this will be getFileName where the hashmaps will be used

        try {
            // Using an HttpResponse class object to catch the response from the server
            HttpResponse response = client.execute(request);
            // The data sent by the server is obtained in this getEntity() function:
            HttpEntity entity = response.getEntity();

            // Extracting the data from the entity object:
            try( InputStream inputStream = entity.getContent();
                 OutputStream outputStream = new FileOutputStream(outputFile) ) {

                // Perform the copy:
                while (true) {

                    // Read from our input:
                    final int bytesRead = inputStream.read(copyBuffer);
                    if (bytesRead == -1) {
                        break;
                    }

                    // Write to our output:
                    outputStream.write(copyBuffer, 0, bytesRead);
                }
            }
            catch (final IOException ex) {
                throw new UserException("Could not copy file: " + dsURL + " -> " + outputFile, ex);
            }
        }
        catch (final IOException ex) {
            throw new UserException("Could not obtain data from "+ dsURL, ex);
        }

        // Extract data sources if requested:       need to change these variables, these will not be right
        if ( extractDataSourcesAfterDownload ) {
            final Path outputDestination = getOutputLocation(dsPath);
            IOUtils.extractTarGz(outputDestination, outputDestination.getParent(), overwriteOutputFile);
        }
        else {
            logger.info("IMPORTANT: You must unzip the downloaded data sources prior to using them with Funcotator.");
        }
    }


    private Path getOutputLocation(final Path dataSourcesPath) {
        if (outputFile == null) {
            return IOUtils.getPath(dataSourcesPath.getFileName().toString());
        } else {
            return outputFile.toPath();
        }
    }

    //==================================================================================================================
    // Helper Data Types:
}



