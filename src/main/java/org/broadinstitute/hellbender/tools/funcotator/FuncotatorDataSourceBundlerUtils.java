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

    private static  Map<String, String> bacteriaMap = new LinkedHashMap<>();
    private static  Map<String, String> fungiMap = new LinkedHashMap<>();
    private static  Map<String, String> metazoaMap = new LinkedHashMap<>();
    private static  Map<String, String> plantsMap = new LinkedHashMap<>();
    private static  Map<String, String> protistsMap = new LinkedHashMap<>();

    private static final String bacteriaFileName = "uniprot_report_EnsemblBacteria_compressed";
    private static final String fungiFileName = "uniprot_report_EnsemblFungi_compressed";
    private static final String metazoaFileName = "uniprot_report_EnsemblMetazoa_compressed";
    private static final String plantsFileName = "uniprot_report_EnsemblPlants_compressed";
    private static final String protistsFileName = "uniprot_report_EnsemblProtists_compressed";

    public static final String ENSEMBL_VERSION = "51";
    public static final String SEPARATOR_CHARACTER = ".";
    public static final String GTF_GZ_ESTENSION = "gtf.gz";

    public static Map<String, String> orgNamesAndFileNames = new LinkedHashMap<>();
    public static final List<String> orgNameKeys = new ArrayList<>();
    public static final List<String> fileNameValues = new ArrayList<>();
    public static Map<String, Map<String, String>> orgNamesMapNames = new LinkedHashMap<>();
    public static final List<Map<String, String>> mapNameValues = new ArrayList<>();

    /**
     * Return the file name associated with the given species and organism.
     * @param orgName
     * @param speciesName
     */
    public static void getDSFileName(String orgName, String speciesName) {

    }

    /**
     * Initialize our hashmaps of lookup tables:
     */
    public static void buildMaps() {
        orgNameKeys.add("bacteria");
        orgNameKeys.add("fungi");
        orgNameKeys.add("metazoa");
        orgNameKeys.add("plants");
        orgNameKeys.add("protists");

        fileNameValues.add(bacteriaFileName);
        fileNameValues.add(fungiFileName);
        fileNameValues.add(metazoaFileName);
        fileNameValues.add(plantsFileName);
        fileNameValues.add(protistsFileName);

        mapNameValues.add(bacteriaMap);
        mapNameValues.add(fungiMap);
        mapNameValues.add(metazoaMap);
        mapNameValues.add(plantsMap);
        mapNameValues.add(protistsMap);


        orgNamesAndFileNames = FuncotatorUtils.createLinkedHashMapFromLists(orgNameKeys, fileNameValues);
        orgNamesMapNames = FuncotatorUtils.createLinkedHashMapFromLists(orgNameKeys, mapNameValues);

        for (String orgName : orgNameKeys) {
            String urlName =  DataSourceUtils.DATA_SOURCES_BASE_URL + DataSourceUtils.DATA_SOURCES_VERSION + orgName + "/" + orgNamesAndFileNames.get(orgName);
            readUniprotFile(urlName, orgName);
        }
    }

    public static void readUniprotFile(String urlFilePath, String orgName){
        final List<String> keys = new ArrayList<>();
        final List<String> values = new ArrayList<>();
        try {
            Scanner inputStream = new Scanner(new File(urlFilePath));
            while (inputStream.hasNextLine()) {
                String data = inputStream.nextLine();
                String[] columnValues = data.split("\t");
                keys.add(columnValues[0].replaceAll("[[]()]]", ""));
                values.add(columnValues[1] + "." + columnValues[4] + "." + ENSEMBL_VERSION);
            }
            inputStream.close();
            if ( orgNamesMapNames.get(orgName) != null ) {
                orgNamesMapNames.get(orgName) = FuncotatorUtils.createLinkedHashMapFromLists(keys, values);
            }
            else {
                throw new UserException("Key has a null value.");
            }

        }
        catch ( final IOException ex ){
            throw new UserException("Could not open file: " + urlFilePath, ex);
        }

    }
}
