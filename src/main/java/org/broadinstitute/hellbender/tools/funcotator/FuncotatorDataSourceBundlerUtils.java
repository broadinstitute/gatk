package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.nio.file.Path;
import java.util.*;
//import org.apache.poi.hssf.usermodel.HSSFSheet;
//import org.apache.poi.hssf.usermodel.HSSFWorkbook;
//import org.apache.poi.ss.usermodel.Cell;
//import org.apache.poi.ss.usermodel.FormulaEvaluator;
//import org.apache.poi.ss.usermodel.Row;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.File;
import java.io.IOException;

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

    // File names for the files which will be used to
    private static final String bacteriaFileName = "uniprot_report_EnsemblBacteria.txt";
    private static final String fungiFileName = "uniprot_report_EnsemblFungi.txt";
    private static final String metazoaFileName = "uniprot_report_EnsemblMetazoa.txt";
    private static final String plantsFileName = "uniprot_report_EnsemblPlants.txt";
    private static final String protistsFileName = "uniprot_report_EnsemblProtists.txt";

    public static final String ENSEMBL_VERSION = "51";
    public static final String SEPARATOR_CHARACTER = ".";

    public static Map<String, String> orgNamesAndFileNames = new LinkedHashMap<>();
    public static Map<String, Map<String, String>> orgNamesMapNames = new LinkedHashMap<>();
    public static final List<String> orgNameKeys = new ArrayList<>();
    public static final List<String> fileNameValues = new ArrayList<>();
    public static final List<Map<String, String>> mapNameValues = new ArrayList<>();

    /**
     * Return the file name associated with the given species and organism.
     * @param orgName
     * @param speciesName
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
     * Initialize our hashmaps of lookup tables:
     */
    public static void buildMaps() {

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
        mapNameValues.add(bacteriaMap);
        mapNameValues.add(fungiMap);
        mapNameValues.add(metazoaMap);
        mapNameValues.add(plantsMap);
        mapNameValues.add(protistsMap);

        // Creating map from the organism names to their corresponding file names:
        orgNamesAndFileNames = FuncotatorUtils.createLinkedHashMapFromLists(orgNameKeys, fileNameValues);
        // Creating map from the organism names to their corresponding map objects:
        orgNamesMapNames = FuncotatorUtils.createLinkedHashMapFromLists(orgNameKeys, mapNameValues);

        // Looping through each organism in our list of valid organisms to initialize their maps from species name to species file name:
        for (String orgName : orgNameKeys) {
            String urlName =  DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
            readUniprotFile(urlName, orgName);
        }

        haveInitializedMaps = true;
    }

    public static void buildMap(String orgName){
        String urlName =  DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
        readUniprotFile(urlName, orgName);
    }

    public static String buildMapGetFileName(String orgName, String speciesName) {
        String urlName = DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
        readUniprotFile(urlName, orgName);

        // Get the map for this organism:
        Map<String, String> tempMap = orgNamesMapNames.get(orgName);

        // Check to see if specified species name is valid and if so, return the corresponding file name:
        if (tempMap.get(speciesName) == null) {
            throw new UserException.BadInput("Given species name: " + speciesName + " is not a valid species for organism: " + orgName + "!");
        } else {
            return tempMap.get(speciesName);

        }
    }

    public static void readUniprotFile(String urlFilePath, String orgName){

        // Keys will be a list of all of the species for a given organism:
        final List<String> keys = new ArrayList<>();
        // Values will be a list of all of the associated file names for each species:
        final List<String> values = new ArrayList<>();

        // Using Scanner object to read in the data on the web page specified by the given url:
        try {
            Scanner inputStream = new Scanner(new File(urlFilePath));
            while (inputStream.hasNextLine()) {
                String data = inputStream.nextLine();
                String[] columnValues = data.split("\t");
                // First column is the species name:
                keys.add(columnValues[0].replaceAll("[[]()]]", ""));
                // Second column is the species division, and third column is the assembly ID, both of which
                // are used to make the corresponding file name for this species:
                values.add(columnValues[1] + "." + columnValues[4] + "." + ENSEMBL_VERSION);
            }
            inputStream.close();

            // Finding the corresponding map for this organism name:
            Map<String, String> tempMap = orgNamesMapNames.get(orgName);

            // Initializing this organism's map from all of its species to their corresponding file names:
            tempMap = FuncotatorUtils.createLinkedHashMapFromLists(keys, values);

            // Setting the organism's map equal to the temporary map value:
//          orgNamesMapNames.get(orgName) = tempMap;
        }
        catch ( final IOException ex ){
            throw new UserException("Could not open file: " + urlFilePath, ex);
        }

    }
}
