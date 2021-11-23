package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.codecs.ProgressReportingDelegatingCodec;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.net.URL;


/**
 * Utilities for reading / working with / manipulating data sources for other organisms.
 * Created to be used with {@link FuncotatorDataSourceBundler}.
 * Created by Hailey on 8/2/21.
 */
public class FuncotatorDataSourceBundlerUtils {
    private static final Logger logger = LogManager.getLogger(FuncotatorDataSourceBundlerUtils.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String ENSEMBL_VERSION                          = "51";
    public static final String SEPARATOR_CHARACTER                      = ".";
    public static final String SEQUENCE_TYPE_AND_STATUS                 = "cdna.all";

    public static final int NUM_UNIPROT_COLUMNS = 11;

    private static final Map<String, Map<String, String>> orgNamesMapNames     = new LinkedHashMap<>();

    //==================================================================================================================
    // Static initializer:

    // Initialize our maps:
    static {
        initializeListsAndMaps();
    }

    /**
     * @return The current date in numerical string format (YearMonthDay / yyyyMMdd).
     */
    public static String getCurrentDateString(){
        return new SimpleDateFormat("yyyyMMdd").format(Calendar.getInstance().getTime());
    }

    /**
     * Builds a map for the specified organism and then returns the file name associated with the given species name.
     * @param kingdom {@link FuncotatorDataSourceBundler.OrganismKingdom} representing the chosen species' kingdom.
     * @param speciesName The specific species we want to know the file name for.
     * @return The file name associated with the given species name.
     */
    public static String getDatasourceBaseName(final FuncotatorDataSourceBundler.OrganismKingdom kingdom, final String speciesName) {

        final String urlName = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + kingdom.toString()
                + "/" + kingdom.getUniprotFileName();
        // Build map from species names to file names:
        readUniprotFile(urlName, kingdom.toString(), false);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if ( orgNamesMapNames.get(kingdom.toString()).get(speciesName) == null ) {
            throw new UserException.BadInput("Given species name: " + speciesName +
                    " is not a valid species for organism: " + kingdom + "!");
        } else{
            return orgNamesMapNames.get(kingdom.toString()).get(speciesName);
        }
    }

    /**
     * Builds a map for the specified organism and then returns the file name associated with the given species name.
     * @param kingdom {@link FuncotatorDataSourceBundler.OrganismKingdom} representing the chosen species' kingdom.
     * @param speciesName The specific species we want to know the file name for.
     * @return The file name associated with the given species name.
     */
    public static String getFastaFileName(final FuncotatorDataSourceBundler.OrganismKingdom kingdom, final String speciesName) {

        final String urlName = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + kingdom.toString()
                + "/" + kingdom.getUniprotFileName();

        // Build map from species names to file names:
        // NOTE: This is getting the data for the FASTA file here with the third argument:
        readUniprotFile(urlName, kingdom.toString(), true);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if ( orgNamesMapNames.get(kingdom.toString()).get(speciesName) == null ) {
            throw new UserException.BadInput("Given species name: " + speciesName +
                    " is not a valid species for organism: " + kingdom + "!");
        } else{
            return orgNamesMapNames.get(kingdom.toString()).get(speciesName);
        }
    }

    /**
     * Initializes a list of organism names, a list of the uniprot_report file names, and a list of empty maps.
     * Initializes a map from the organism names to their associated uniprot file name.
     * Initializes a map from the organism names to their associated map objects.
     */
    public static void initializeListsAndMaps() {
        for ( final FuncotatorDataSourceBundler.OrganismKingdom k : FuncotatorDataSourceBundler.OrganismKingdom.values() ) {
            orgNamesMapNames.put(k.toString(), new LinkedHashMap<>());
        }
    }

    /**
     * Scans uniprot report file to build file names for each species in the file.
     * @param urlFilePath Url at which the uniprot report file is found.
     * @param orgName Name of the chosen organism.
     * @param isFasta Boolean indicating if we are downloading the fasta file.
     */
    public static void readUniprotFile(final String urlFilePath, final String orgName, final boolean isFasta){

        // Keys will be a list of all of the species for a given organism:
        final List<String> keys = new ArrayList<>();
        // Values will be a list of all of the associated file names for each species:
        final List<String> values = new ArrayList<>();
        // Filename variable
        String fileName;

        try {
            URL uniprotURL = new URL(urlFilePath);
            // Using Scanner object to read in the data on the web page specified by the given url:
            try (final Scanner inputStream = new Scanner(uniprotURL.openStream())){
                while (inputStream.hasNextLine()) {
                    String data = inputStream.nextLine();
                    String[] columnValues = data.split("\t");

                    // Make sure we have the right number of columns:
                    if (columnValues.length != NUM_UNIPROT_COLUMNS ) {
                        throw new GATKException("Uniprot file contains unexpected number of columns: " + columnValues.length + " != " + NUM_UNIPROT_COLUMNS);
                    }

                    // Second column is the species name, these will be the keys in the map from species name to file name:
                    keys.add(columnValues[1]);

                    // If first character is _, then the first letter is not capitalized:
                    if ( columnValues[1].charAt(0) == '_' ) {
                        // Sixth column is the assembly build name, so format is <species name>.<assembly name>
                        fileName = columnValues[1] + "." + columnValues[5];
                    }
                    else {
                        // Otherwise, we have to capitalize first letter of file name:
                        fileName = columnValues[1].substring(0, 1).toUpperCase() + columnValues[1].substring(1) + "." + columnValues[5];
                    }
                    // If we are building the map for fasta files, we have to add the sequence type and status to the file name:
                    if (isFasta) {
                        values.add(fileName + "." + SEQUENCE_TYPE_AND_STATUS);
                    } else {
                        // Otherwise, add the ensembl version number for regular data files:
                        values.add(fileName + "." + ENSEMBL_VERSION);
                    }
                }

                // Creating new linked hash map using keys and values lists and then putting this map as the value for the
                // correct organism:
                orgNamesMapNames.put(orgName, FuncotatorUtils.createLinkedHashMapFromLists(keys, values));
            }
            catch ( final IOException ex ){
                throw new UserException("Could not open file: " + urlFilePath, ex);
            }
        } catch (MalformedURLException ex) {
            throw new UserException("Unable to access URL " + urlFilePath + "!", ex);
        }
    }

    /**
     * Extracts the .gz file given by {@code gzFilePath}.
     * Input {@link String} MUST be to a gzipped gtf or fasta file.
     * Will extract contents in the containing folder of {@code gzFilePath}.
     * Will throw an exception if files exist already.
     * @param gzFilePath {@link String} path to a gzipped gtf or fasta file for extraction.
     * @param decompressedFilePath {@link String} path where the unzipped gtf or fasta file will go.
     * @param doOverwrite {@link boolean} which determine if a pre-existing file should be overwritten or not.
     */
    public static void extractGzFile(final String gzFilePath, final String decompressedFilePath, final boolean doOverwrite) {
        IOUtils.ensurePathIsOkForOutput(IOUtils.getPath(decompressedFilePath), doOverwrite);
        IOUtils.gunzip(IOUtils.getPath(gzFilePath).toFile(), IOUtils.getPath(decompressedFilePath).toFile());
    }

    /**
     * Create an index for a GTF file.
     * @param featurePath {@link Path} of the GTF file to be indexed.
     * @param idxFilePath {@link Path} to which to write the index.
     */
    public static void indexGTF(final Path featurePath, final Path idxFilePath) {
        if (!Files.isReadable(featurePath)) {
            throw new UserException.CouldNotReadInputFile(featurePath);
        }
        // Get the right codec for the file to be indexed. This call will throw an appropriate exception
        // if featureFile is not in a supported format or is unreadable.
        final FeatureCodec<? extends Feature, ?> codec = new ProgressReportingDelegatingCodec<>(
                FeatureManager.getCodecForFile(featurePath), ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES);
        try {
            final Index index = IndexFactory.createDynamicIndex(featurePath, codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
            logger.info("Creating Tribble index.");
            final Path indexPath = Tribble.indexPath(featurePath);
            try {
                index.write(indexPath);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile("Could not write index to file " + indexPath, e);
            }
        } catch ( TribbleException e) {
            // Underlying cause here is usually a malformed file, but can also be things like
            // "codec does not support tabix"
            throw new UserException.CouldNotIndexFile(featurePath, e);
        }
        logger.info("Successfully wrote index to " + idxFilePath);
    }

}
