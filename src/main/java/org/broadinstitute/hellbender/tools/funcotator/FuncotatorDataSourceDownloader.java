package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.codec.digest.MessageDigestAlgorithms;
import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.ArchiveInputStream;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeter;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeterResults;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.nio.file.Files;
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
public class FuncotatorDataSourceDownloader extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(FuncotatorDataSourceDownloader.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String VALIDATE_INTEGRITY_ARG_LONG_NAME = "validate-integrity";
    public static final String SOMATIC_ARG_LONG_NAME            = "somatic";
    public static final String GERMLINE_ARG_LONG_NAME           = "germline";
    public static final String OVERWRITE_ARG_LONG_NAME          = "overwrite-output-file";
    public static final String EXTRACT_AFTER_DOWNLOAD           = "extract-after-download";

    //==================================================================================================================
    // Private Static Members:

    private static String BASE_URL = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.4.20180829";

    private static String GERMLINE_GCLOUD_DATASOURCES_BASEURL     = BASE_URL + "g";
    private static Path   GERMLINE_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(GERMLINE_GCLOUD_DATASOURCES_BASEURL + ".tar.gz");
    private static Path   GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(GERMLINE_GCLOUD_DATASOURCES_BASEURL + ".sha256");

    private static String SOMATIC_GCLOUD_DATASOURCES_BASEURL     = BASE_URL + "s";
    private static Path   SOMATIC_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(SOMATIC_GCLOUD_DATASOURCES_BASEURL + ".tar.gz");
    private static Path   SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(SOMATIC_GCLOUD_DATASOURCES_BASEURL + ".sha256");

    //==================================================================================================================
    // Private Members:

    @Argument(fullName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            shortName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            doc = "Validate the integrity of the data sources after downloading them using sha256.",
            optional = true)
    private boolean doValidateIntegrity = false;

    @Argument(fullName = SOMATIC_ARG_LONG_NAME,
            shortName = SOMATIC_ARG_LONG_NAME,
            mutex = {GERMLINE_ARG_LONG_NAME},
            doc = "Download the latest pre-packaged datasources for somatic functional annotation.",
            optional = true)
    private boolean getSomaticDataSources = false;

    @Argument(fullName = GERMLINE_ARG_LONG_NAME,
            shortName = GERMLINE_ARG_LONG_NAME,
            mutex = {SOMATIC_ARG_LONG_NAME},
            doc = "Download the latest pre-packaged datasources for germline functional annotation.",
            optional = true)
    private boolean getGermlineDataSources = false;

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName = OVERWRITE_ARG_LONG_NAME,
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

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

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

    @Override
    protected Object doWork() {

        // Get the correct data source:
        if ( getSomaticDataSources ) {
            downloadAndValidateDataSources("Somatic", SOMATIC_GCLOUD_DATASOURCES_PATH, SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH);
        }

        if ( getGermlineDataSources ) {
            downloadAndValidateDataSources("Germline", GERMLINE_GCLOUD_DATASOURCES_PATH, GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH);
        }

        // Token return value:
        return true;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private void extractDataSources(final Path localDataSourcesPath) {

        logger.info("Extracting data sources archive: " + localDataSourcesPath.toUri());

        // Create a stream for the data sources input.
        // (We know it will be a tar.gz):
        try ( final InputStream fi = Files.newInputStream(localDataSourcesPath);
              final InputStream bi = new BufferedInputStream(fi);
              final InputStream gzi = new GzipCompressorInputStream(bi);
              final ArchiveInputStream archiveStream = new TarArchiveInputStream(gzi)) {

            extractFilesFromArchiveStream(archiveStream, localDataSourcesPath);
        }
        catch (final IOException ex) {
            throw new UserException("Could not extract data from: " + localDataSourcesPath.toUri());
        }

    }

    private void extractFilesFromArchiveStream(final ArchiveInputStream archiveStream,
                                               final Path localDataSourcesPath) throws IOException {

        // Adapted from: http://commons.apache.org/proper/commons-compress/examples.html

        // Go through the archive and get the entries:
        ArchiveEntry entry;
        while ((entry = archiveStream.getNextEntry()) != null) {

            logger.info("Extracting file: " + entry.getName());

            // Make sure we can read the data for the entry:
            if (!archiveStream.canReadEntryData(entry)) {
                throw new UserException("Could not read data from archive file(" + localDataSourcesPath.toUri() + "): " + entry.getName());
            }

            // Get the path for the entry on disk and make sure it's OK:
            final Path extractedEntryPath = localDataSourcesPath.resolveSibling(IOUtils.getPath(entry.getName()));
            ensurePathIsOkForOutput(extractedEntryPath);

            // Now we can create the entry in our output location:
            final File extractedEntryFile = extractedEntryPath.toFile();
            if (entry.isDirectory()) {
                // Handle a directory entry:
                if (!extractedEntryFile.isDirectory() && !extractedEntryFile.mkdirs()) {
                    throw new IOException("Failed to create directory " + extractedEntryFile);
                }
            } else {
                // Handle a file entry:

                // Make sure we can actually create the file in the expected/requested directory:
                final File parent = extractedEntryFile.getParentFile();
                if (!parent.isDirectory() && !parent.mkdirs()) {
                    throw new IOException("Failed to create directory " + parent);
                }

                // Create the output file from the stream:
                try (final OutputStream o = Files.newOutputStream(extractedEntryPath)) {
                    org.apache.commons.io.IOUtils.copy(archiveStream, o);
                }
            }
        }
    }

    private void ensurePathIsOkForOutput(final Path p) {
        if ( Files.exists(p) ) {
            if ( overwriteOutputFile ) {
                logger.warn("Overwriting existing output destination: " + p.toUri());
            }
            else {
                throw new UserException("Output destination already exists: " + p.toUri());
            }
        }
    }

    private NioFileCopierWithProgressMeter createNioDownloader(final Path dsPath) {
        // Get the data sources file:
        final Path outputDestination = getOutputLocation(dsPath);

        // Set up and initiate our download:
        return NioFileCopierWithProgressMeter.create(dsPath, outputDestination, overwriteOutputFile);
    }

    private void downloadAndValidateDataSources(final String dataSourceType, final Path dsPath, final Path checksumPath) {
        logger.info(dataSourceType + " data sources selected.");

        // Create our downloader:
        final NioFileCopierWithProgressMeter xerox = createNioDownloader(dsPath);

        // Get the datasources file:
        final NioFileCopierWithProgressMeterResults results = downloadDataSources(xerox, checksumPath);

        // Confirm file integrity if requested:
        if ( doValidateIntegrity ) {
            validateIntegrity(results);
        }

        // Extract data sources if requested:
        if ( extractDataSourcesAfterDownload ) {
            extractDataSources(results.getDestination());
        }
        else {
            logger.info("IMPORTANT: You must unzip the downloaded data sources prior to using them with Funcotator.");
        }
    }

    private NioFileCopierWithProgressMeterResults downloadDataSources(final NioFileCopierWithProgressMeter xerox,
                                                                      final Path checksumPath) {

        // Set up the validity check while the download occurs if requested:
        if ( doValidateIntegrity ) {
            // Read the sha256sum into memory:
            final String expectedSha256Sum = readSha256SumFromPath(checksumPath);

            // Setup the copier to calculate the checksum:
            xerox.setChecksumAlgorithmAndExpectedChecksum(MessageDigestAlgorithms.SHA_256, expectedSha256Sum);
        }

        // Initiate the copy and return our results:
        return xerox.initiateCopy();
    }

    private void validateIntegrity(final NioFileCopierWithProgressMeterResults results) {

        // verify the hashes are the same:
        if ( !results.isDestFileValid() ) {
            throw new UserException("ERROR: downloaded data sources are corrupt!  Unexpected checksum: " + results.getChecksum() + " != " + results.getExpectedChecksum());
        }
        else {
            logger.info("Integrity check on downloaded data sources succeeded.");
        }
    }

    private String readSha256SumFromPath(final Path sha256Path) {
        final String expectedSha256Sum;
        try {
            logger.info("Collecting expected checksum...");
            expectedSha256Sum = Files.lines(sha256Path).findFirst().orElse(null);

            if ( expectedSha256Sum == null ) {
                throw new UserException("Unable to retrieve expected checksum from: " + sha256Path.toUri());
            }

            logger.info("Collection complete!");
        }
        catch ( final IOException ex ) {
            throw new UserException("Could not read in sha256sum from file: " + sha256Path.toUri().toString(), ex);
        }

        // Clean up and return the checksum:
        return cleanExpectedSha256SumString(expectedSha256Sum);
    }

    private Path getOutputLocation(final Path dataSourcesPath) {
        if ( outputFile == null ) {
            return IOUtils.getPath(dataSourcesPath.getFileName().toString());
        }
        else {
            return outputFile.toPath();
        }
    }

    private String cleanExpectedSha256SumString(final String expectedSha256SumString) {

        String cleanString = expectedSha256SumString.trim().toLowerCase();

        // The format of the file can contain the checksum, followed by the file name.
        // If this is the case, we need to truncate the string:
        if ( cleanString.contains(" ") ) {
            cleanString = cleanString.substring(0, cleanString.indexOf(" "));
        }
        if ( cleanString.contains("\t")) {
            cleanString = cleanString.substring(0, cleanString.indexOf("\t"));
        }

        return cleanString;
    }

    //==================================================================================================================
    // Helper Data Types:

}
