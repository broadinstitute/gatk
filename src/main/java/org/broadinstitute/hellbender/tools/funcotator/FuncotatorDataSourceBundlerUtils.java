package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;

import java.net.MalformedURLException;
import java.util.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import org.broadinstitute.hellbender.utils.io.IOUtils;


/**
 * Utilities for reading / working with / manipulating data sources for other organisms.
 * Created to be used with {@link FuncotatorDataSourceBundler}.
 * Created by Hailey
 */
public class FuncotatorDataSourceBundlerUtils {

    private final static Logger logger = LogManager.getLogger(FuncotatorDataSourceBundlerUtils.class);

    // If we have already initialized the maps then there is no need to build them again:
    public static boolean haveInitializedMaps = false;

    // Maps will map from each species in a given organism to the species file name associated with that species name:
    private static  Map<String, String> bacteriaMap = new LinkedHashMap<>();
    private static  Map<String, String> fungiMap = new LinkedHashMap<>();
    private static  Map<String, String> metazoaMap = new LinkedHashMap<>();
    private static  Map<String, String> plantsMap = new LinkedHashMap<>();
    private static  Map<String, String> protistsMap = new LinkedHashMap<>();

    // File names for the files which will be used to initialize the maps from species name to file name
    private static final String bacteriaFileName = "uniprot_report_EnsemblBacteria.txt";
    private static final String fungiFileName = "uniprot_report_EnsemblFungi.txt";
    private static final String metazoaFileName = "uniprot_report_EnsemblMetazoa.txt";
    private static final String plantsFileName = "uniprot_report_EnsemblPlants.txt";
    private static final String protistsFileName = "uniprot_report_EnsemblProtists.txt";

    public static final String ENSEMBL_VERSION = "51";
    public static final String SEPARATOR_CHARACTER = ".";
    public static final String SEQUENCE_TYPE_AND_STATUS = "cdna.all";

    public static Map<String, String> orgNamesAndFileNames = new LinkedHashMap<>();
    public static Map<String, Map<String, String>> orgNamesMapNames = new LinkedHashMap<>();
    public static final List<String> orgNameKeys = new ArrayList<>();
    public static final List<String> fileNameValues = new ArrayList<>();
    public static final List<Map<String, String>> nullMapValues = new ArrayList<>();

    /**
     * Return the file name associated with the given species and organism.
     * Note: This function uses the buildMaps helper function which builds all five maps
     * instead of one, so this function takes longer than buildMapGetFileName.
     * @param orgName Name of the chosen organism.
     * @param speciesName Name of the chosen species.
     */
    public static String getDSFileName(String orgName, String speciesName) {

        // Initialize the maps from species name to species file name for each organism:
        buildMaps();

        // Get the map for this organism:
        Map<String, String> tempMap = orgNamesMapNames.get(orgName);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if (tempMap.get(speciesName) == null) {
            throw new UserException.BadInput("Given species name: " + speciesName + " is not a valid species for organism: " + orgName + "!");
        }
        else {
            return tempMap.get(speciesName);
        }
    }

    /**
     * Builds all five maps from species name to file name for the regular data sources.
     */
    public static void buildMaps() {

        buildListsAndMaps();

        // Looping through each organism in our list of valid organisms to initialize their maps from species name to species file name:
        for (String orgName : orgNameKeys) {
            String urlName =  DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
            readUniprotFile(urlName, orgName, false);
        }

        haveInitializedMaps = true;
    }


    /**
     * Builds a map for the specified organism and then returns the file name associated with the given species name.
     * @param orgName The chosen organism.
     * @param speciesName The specific species we want to know the file name for.
     * @param isFasta Boolean to determine if we are building a Fasta map or a regular map.
     * @return
     */
    public static String buildMapGetFileName(String orgName, String speciesName, boolean isFasta) {

        // Only need to initialize these variables if we have not done it yet
        if (!isFasta) {
            buildListsAndMaps();
        }

        String urlName = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
        // Build map from species names to file names:
        readUniprotFile(urlName, orgName, isFasta);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if ( orgNamesMapNames.get(orgName).get(speciesName) == null ) {
            throw new UserException.BadInput("Given species name: " + speciesName + " is not a valid species for organism: " + orgName + "!");
        } else{
            return orgNamesMapNames.get(orgName).get(speciesName);
        }
    }

    /**
     * Initializes a list of organism names, a list of the uniprot_report file names, and a list of empty maps.
     * Initializes a map from the organism names to their associated uniprot file name.
     * Initializes a map from the organism names to their associated map objects.
     */
    public static void buildListsAndMaps() {
        // Initializing list of allowed organism names:
        orgNameKeys.add("bacteria");
        orgNameKeys.add("fungi");
        orgNameKeys.add("metazoa");
        orgNameKeys.add("plants");
        orgNameKeys.add("protists");

        // Initializing list of file names for each organism:
        fileNameValues.add(bacteriaFileName);
        fileNameValues.add(fungiFileName);
        fileNameValues.add(metazoaFileName);
        fileNameValues.add(plantsFileName);
        fileNameValues.add(protistsFileName);

        // Initializing list of map names for each organism:
        nullMapValues.add(bacteriaMap);
        nullMapValues.add(fungiMap);
        nullMapValues.add(metazoaMap);
        nullMapValues.add(plantsMap);
        nullMapValues.add(protistsMap);

        // Creating map from the organism names to their corresponding file names:
        orgNamesAndFileNames = FuncotatorUtils.createLinkedHashMapFromLists(orgNameKeys, fileNameValues);
        // Creating map from the organism names to their corresponding map objects:
        orgNamesMapNames = FuncotatorUtils.createLinkedHashMapFromLists(orgNameKeys, nullMapValues);
    }

    /**
     * Scans uniprot report file to build file names for each species in the file.
     * @param urlFilePath Url at which the uniprot report file is found.
     * @param orgName Name of the chosen organism.
     * @param isFasta Boolean indicating if we are downloading the fasta file.
     */
    public static void readUniprotFile(String urlFilePath, String orgName, boolean isFasta){

        // Keys will be a list of all of the species for a given organism:
        final List<String> keys = new ArrayList<>();
        // Values will be a list of all of the associated file names for each species:
        final List<String> values = new ArrayList<>();
        // Filename variable
        String fileName;

        try {
            URL uniprotURL = new URL(urlFilePath);
            // Using Scanner object to read in the data on the web page specified by the given url:
            try {
                Scanner inputStream = new Scanner(uniprotURL.openStream());
                while (inputStream.hasNextLine()) {
                    String data = inputStream.nextLine();
                    String[] columnValues = data.split("\t");
                    // Second column is the species name, these will be the keys in the map from species name to file name:
                    keys.add(columnValues[1]);

                    // If first character is _, then the first letter is not capitalized:
                    if (columnValues[1].substring(0, 1).equals("_")) {
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
                inputStream.close();

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
     * Extracts the gtf.gz file given by {@code gtfGzFilePath}.
     * Input {@link String} MUST be to a gzipped gtf file.
     * Will extract contents in the containing folder of {@code gtfGzFilePath}.
     * Will throw an exception if files exist already.
     * @param gtfGzFilePath {@link String} to a gzipped gtf file for extraction.
     */
    public static void extractGtfGz(String gtfGzFilePath, String decompressedFilePath, boolean doOverwrite) {
        byte[] buffer = new byte[1024];

        IOUtils.ensurePathIsOkForOutput(IOUtils.getPath(decompressedFilePath), doOverwrite);
        try {
            FileInputStream inputStream = new FileInputStream(gtfGzFilePath);
            GZIPInputStream gZIPInputStream = new GZIPInputStream(inputStream);
            FileOutputStream fileOutputStream = new FileOutputStream(decompressedFilePath);

            int bytes_read;

            while ((bytes_read = gZIPInputStream.read(buffer)) > 0) {
                fileOutputStream.write(buffer, 0, bytes_read);
            }
        } catch (IOException ex) {
            throw new UserException("Could not obtain data from " + gtfGzFilePath, ex);
        }

    }
}
