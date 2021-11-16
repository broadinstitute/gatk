package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.net.MalformedURLException;
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

    //==================================================================================================================
    // Private Static Members:

//    // If we have already initialized the maps then there is no need to build them again:
//    public static boolean haveInitializedMaps = false;

    // These map each species for a given organism to the file name associated with that species name:
    private static final Map<String, String> bacteriaMap     = new LinkedHashMap<>();
    private static final Map<String, String> fungiMap        = new LinkedHashMap<>();
    private static final Map<String, String> metazoaMap      = new LinkedHashMap<>();
    private static final Map<String, String> plantsMap       = new LinkedHashMap<>();
    private static final Map<String, String> protistsMap     = new LinkedHashMap<>();

    // File names for the files which will be used to initialize the maps from species name to file name:
    private static final String bacteriaFileName        = "uniprot_report_EnsemblBacteria.txt";
    private static final String fungiFileName           = "uniprot_report_EnsemblFungi.txt";
    private static final String metazoaFileName         = "uniprot_report_EnsemblMetazoa.txt";
    private static final String plantsFileName          = "uniprot_report_EnsemblPlants.txt";
    private static final String protistsFileName        = "uniprot_report_EnsemblProtists.txt";

    //==================================================================================================================
    // Public Static Members:

    public static final String ENSEMBL_VERSION                          = "51";
    public static final String SEPARATOR_CHARACTER                      = ".";
    public static final String SEQUENCE_TYPE_AND_STATUS                 = "cdna.all";

    public static Map<String, String> orgNamesAndFileNames              = new LinkedHashMap<>();
    public static Map<String, Map<String, String>> orgNamesMapNames     = new LinkedHashMap<>();

    public static final List<String> orgNameKeys                        = new ArrayList<>();
    public static final List<String> fileNameValues                     = new ArrayList<>();
    public static final List<Map<String, String>> nullMapValues         = new ArrayList<>();

    //==================================================================================================================
    // Static initializer:

    // Initialize our maps:
    static {
        initializeListsAndMaps();
    }

    //==================================================================================================================
    // Public Static Methods:

//    /**
//     * Return the file name associated with the given species and organism.
//     * Note: This function uses the buildMaps helper function which builds all five maps
//     * instead of one, so this function takes longer than buildMapGetFileName.
//     * @param orgName Name of the chosen organism.
//     * @param speciesName Name of the chosen species.
//     */
//    public static String getDSFileName(final String orgName, final String speciesName) {
//
//        // Initialize the maps from species name to species file name for each organism:
//        buildMaps();
//
//        // Get the map for this organism:
//        Map<String, String> tempMap = orgNamesMapNames.get(orgName);
//
//        // Check to see if specified species name is valid and if so, return the corresponding file name:
//        if (tempMap.get(speciesName) == null) {
//            throw new UserException.BadInput("Given species name: " + speciesName + " is not a valid species for organism: " + orgName + "!");
//        }
//        else {
//            return tempMap.get(speciesName);
//        }
//    }

//    /**
//     * Builds all five maps from species name to file name for the regular data sources.
//     */
//    public static void buildMaps() {
//
//        buildListsAndMaps();
//
//        // Looping through each organism in our list of valid organisms to initialize their maps from species name to species file name:
//        for (final String orgName : orgNameKeys) {
//            String urlName =  DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
//            readUniprotFile(urlName, orgName, false);
//        }
//
//        haveInitializedMaps = true;
//    }

    /**
     * @return The current date in numerical string format (YearMonthDay / yyyyMMdd).
     */
    public static String getCurrentDateString(){
        return new SimpleDateFormat("yyyyMMdd").format(Calendar.getInstance().getTime());
    }

    /**
     * Builds a map for the specified organism and then returns the file name associated with the given species name.
     * @param orgName The chosen organism.
     * @param speciesName The specific species we want to know the file name for.
     * @return The file name associated with the given species name.
     */
    public static String getDatasourceBaseName(final String orgName, final String speciesName) {

        final String urlName = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName
                + "/" + orgNamesAndFileNames.get(orgName);
        // Build map from species names to file names:
        readUniprotFile(urlName, orgName, false);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if ( orgNamesMapNames.get(orgName).get(speciesName) == null ) {
            throw new UserException.BadInput("Given species name: " + speciesName +
                    " is not a valid species for organism: " + orgName + "!");
        } else{
            return orgNamesMapNames.get(orgName).get(speciesName);
        }
    }

    /**
     * Builds a map for the specified organism and then returns the file name associated with the given species name.
     * @param orgName The chosen organism.
     * @param speciesName The specific species we want to know the file name for.
     * @return The file name associated with the given species name.
     */
    public static String getFastaFileName(final String orgName, final String speciesName) {

        final String urlName = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName
                + "/" + orgNamesAndFileNames.get(orgName);
        // Build map from species names to file names:
        readUniprotFile(urlName, orgName, true);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if ( orgNamesMapNames.get(orgName).get(speciesName) == null ) {
            throw new UserException.BadInput("Given species name: " + speciesName +
                    " is not a valid species for organism: " + orgName + "!");
        } else{
            return orgNamesMapNames.get(orgName).get(speciesName);
        }
    }

    /**
     * Initializes a list of organism names, a list of the uniprot_report file names, and a list of empty maps.
     * Initializes a map from the organism names to their associated uniprot file name.
     * Initializes a map from the organism names to their associated map objects.
     */
    public static void initializeListsAndMaps() {
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
            try {
                Scanner inputStream = new Scanner(uniprotURL.openStream());
                while (inputStream.hasNextLine()) {
                    String data = inputStream.nextLine();
                    String[] columnValues = data.split("\t");
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

}
