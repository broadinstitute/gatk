package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.NioFileCopier;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Path;

/**
 * Tool to download the latest data sources for {@link Funcotator}.
 * Created by jonn on 8/23/18.
 */
@CommandLineProgramProperties(
        summary = "Download the Broad Institute pre-packaged data sources for the somatic or germline use case for Funcotator.",
        oneLineSummary = "Data source downloader for Funcotator.",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FuncotatorDataSourceRetriever extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(FuncotatorDataSourceRetriever.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String VALIDATE_INTEGRITY_ARG_LONG_NAME = "validate-integrity";
    public static final String SOMATIC_ARG_LONG_NAME            = "somatic";
    public static final String GERMLINE_ARG_LONG_NAME           = "germline";
    public static final String OVERWRITE_ARG_LONG_NAME          = "overwrite-output-file";

    //==================================================================================================================
    // Private Static Members:

    private static String GERMLINE_GCLOUD_DATASOURCES_BASEURL     = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.4.20180829g";
    private static Path   GERMLINE_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(GERMLINE_GCLOUD_DATASOURCES_BASEURL + ".tar.gz");
    private static Path   GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(GERMLINE_GCLOUD_DATASOURCES_BASEURL + ".sha256");

    private static String SOMATIC_GCLOUD_DATASOURCES_BASEURL     = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.4.20180829s";
    private static Path   SOMATIC_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(SOMATIC_GCLOUD_DATASOURCES_BASEURL + ".tar.gz");
    private static Path   SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(SOMATIC_GCLOUD_DATASOURCES_BASEURL + ".sha256");

    //==================================================================================================================
    // Private Members:

    @Argument(fullName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            shortName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            doc = "Validate the integrity of the data sources after downloading them using sha256.", optional = true)
    private boolean validateIntegrity = false;

    @Argument(fullName = SOMATIC_ARG_LONG_NAME,
            shortName = SOMATIC_ARG_LONG_NAME,
            mutex = {GERMLINE_ARG_LONG_NAME},
            doc = "Download the latest pre-packaged datasources for somatic functional annotation.", optional = true)
    private boolean getSomaticDataSources = false;

    @Argument(fullName = GERMLINE_ARG_LONG_NAME,
            shortName = GERMLINE_ARG_LONG_NAME,
            mutex = {SOMATIC_ARG_LONG_NAME},
            doc = "Download the latest pre-packaged datasources for germline functional annotation.", optional = true)
    private boolean getGermlineDataSources = false;

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName = OVERWRITE_ARG_LONG_NAME,
            doc = "Overwrite output file if it exists already.", optional = true)
    private boolean overwriteOutputFile = false;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    protected Object doWork() {

        // Get the correct data source:
        if ( getSomaticDataSources ) {
            logger.info("Somatic data sources selected.");

            // Get the germline datasources file:
            final NioFileCopier xerox = downloadDatasources(SOMATIC_GCLOUD_DATASOURCES_PATH);

            // Validate the data sources if requested:
            if (validateIntegrity) {
                logger.info("Integrity validation selected.");
                validateIntegrity(xerox, SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH);
            }
        }

        if ( getGermlineDataSources ) {

            logger.info("Germline data sources selected.");

            // Get the germline datasources file:
            final NioFileCopier xerox = downloadDatasources(GERMLINE_GCLOUD_DATASOURCES_PATH);

            // Validate the data sources if requested:
            if (validateIntegrity) {
                logger.info("Integrity validation selected.");
                validateIntegrity(xerox, GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH);
            }
        }

        logger.info("IMPORTANT: You must unzip the downloaded data sources prior to using them with Funcotator.");

        // Token return value:
        return true;
    }

    @Override
    protected void onStartup() {
        // Make sure the user specified at least one data source to download:
        if ((!getSomaticDataSources) && (!getGermlineDataSources)) {
            throw new UserException("Must select either somatic or germline datasources.");
        }

        if ( overwriteOutputFile ) {
            logger.info("Overwrite ENABLED.  Will overwrite existing data sources download.");
        }
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private NioFileCopier downloadDatasources(final Path dataSourcesPath) {

        // Get the data sources file:
        final Path outputDestination = getDataSourceLocalPath(dataSourcesPath);

        // Set up and initiate our download:
        final NioFileCopier xerox = NioFileCopier.create(dataSourcesPath, outputDestination, overwriteOutputFile);
        xerox.initiateCopy();

        return xerox;
    }

    private void validateIntegrity(final NioFileCopier fileCopier, final Path remoteSha256Path) {

        // Get the SHA 256 file:
        final Path localSha256Path = getSha256File(remoteSha256Path);

        // Read the sha256sum into memory:
        final String expectedSha256Sum = readSha256SumFromFile(localSha256Path);

        // Calculate the sha256sum of the data sources file:
        final boolean downloadedDataSourcesAreValid = fileCopier.validateIntegrity(expectedSha256Sum, DigestUtils::sha256Hex);

        // verify the hashes are the same:
        if ( !downloadedDataSourcesAreValid ) {
            throw new GATKException("ERROR: downloaded data sources are corrupt!  Unexpected checksum: " + fileCopier.getLatestChecksum() + " != " + expectedSha256Sum);
        }
        else {
            logger.info("Data sources are valid.");
        }
    }

    private Path getSha256File(final Path remoteSha256Path) {
        final Path localSha256Path = getDataSourceLocalPath(remoteSha256Path);
        logger.info("Retrieving expected checksum file...");
        NioFileCopier.create(remoteSha256Path, localSha256Path, overwriteOutputFile).setSilentCopy(true).initiateCopy();
        logger.info("File transfer complete!");

        return localSha256Path;
    }

    private String readSha256SumFromFile(final Path localSha256Path) {
        final String expectedSha256Sum;
        try {
            logger.info("Collecting expected checksum...");
            expectedSha256Sum = FileUtils.readFileToString(localSha256Path.toFile(), Charset.defaultCharset());
            logger.info("Collection complete!");
        }
        catch ( final IOException ex ) {
            throw new GATKException("Could not read in sha256sum from file: " + localSha256Path.toUri().toString(), ex);
        }

        // Clean up and return the checksum:
        return cleanExpectedSha256SumString(expectedSha256Sum);
    }

    private Path getDataSourceLocalPath(final Path dataSourcesPath) {
        return IOUtils.getPath(dataSourcesPath.getFileName().toString());
    }

    private String cleanExpectedSha256SumString(final String expectedSha256SumString) {

        String cleanString = expectedSha256SumString.trim().toLowerCase();

        // The format of the file can contain the checksum, followed by the file name.
        // If this is the case, we need to truncate the string:
        while ( cleanString.contains(" ") ) {
            cleanString = cleanString.substring(0, cleanString.indexOf(" "));
        }
        while (  cleanString.contains("\t") ) {
            cleanString = cleanString.substring(0, cleanString.indexOf("\t"));
        }

        return cleanString;
    }

    //==================================================================================================================
    // Helper Data Types:

}
