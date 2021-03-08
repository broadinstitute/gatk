package org.broadinstitute.hellbender.tools.variantdb;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

public class SampleList {
    static final Logger logger = LogManager.getLogger(SampleList.class);

    private Map<Long, String> sampleIdMap = new HashMap<>();
    private Map<String, Long> sampleNameMap = new HashMap<>();

    public SampleList(String sampleTableName, File sampleFile, boolean printDebugInformation) {
        if (sampleTableName != null) {
            initializeMaps(new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS), printDebugInformation);
        } else if (sampleFile != null) {
            initializeMaps(sampleFile);
        } else {
            throw new IllegalArgumentException("--sample-file or --sample-table must be provided.");
        }
    }

    public int size() {
        return sampleIdMap.size();
    }

    public Collection<String> getSampleNames() {
        return sampleIdMap.values();
    }

    public String getSampleName(long id) {
        return sampleIdMap.get(id);
    }

    public long getSampleId(String name) {
        return sampleNameMap.get(name);
    }

    public Map<Long, String> getMap() {
        return sampleIdMap;
    }

//    protected Map<String, Integer> getSampleNameMap(TableReference sampleTable, List<String> samples, boolean printDebugInformation) {
//        Map<String, Integer> results = new HashMap<>();
//        // create optional where clause
//        String whereClause = "";
//        if (samples != null && samples.size() > 0) {
//            whereClause = " WHERE " + SchemaUtils.SAMPLE_NAME_FIELD_NAME + " in (\'" + StringUtils.join(samples, "\',\'") + "\') ";
//        }
//
//        TableResult queryResults = querySampleTable(sampleTable.getFQTableName(), whereClause, printDebugInformation);
//
//        // Add our samples to our map:
//        for (final FieldValueList row : queryResults.iterateAll()) {
//            results.put(row.get(1).getStringValue(), (int) row.get(0).getLongValue());
//        }
//        return results;
//    }


    protected void initializeMaps(TableReference sampleTable, boolean printDebugInformation) {
        TableResult queryResults = querySampleTable(sampleTable.getFQTableName(), "", printDebugInformation);

        // Add our samples to our map:
        for (final FieldValueList row : queryResults.iterateAll()) {
            long id = row.get(0).getLongValue();
            String name = row.get(1).getStringValue();
            sampleIdMap.put(id, name);
            sampleNameMap.put(name, id);
        }
    }

    protected void initializeMaps(File cohortSampleFile) {
        try {
            Files.readAllLines(cohortSampleFile.toPath(), StandardCharsets.US_ASCII).stream()
                .map(s -> s.split(","))
                .forEach(tokens -> {
                    long id = Long.parseLong(tokens[0]);
                    String name = tokens[1];
                    sampleIdMap.put(id, name);
                    sampleNameMap.put(name, id);
                });
        } catch (IOException e) {
            throw new IllegalArgumentException("Could not parse --cohort-sample-file", e);
        }
    }

    private TableResult querySampleTable(String fqSampleTableName, String whereClause, boolean printDebugInformation) {
        // Get the query string:
        final String sampleListQueryString =
                "SELECT " + SchemaUtils.SAMPLE_ID_FIELD_NAME + ", " + SchemaUtils.SAMPLE_NAME_FIELD_NAME +
                        " FROM `" + fqSampleTableName + "`" + whereClause;


        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        // Show our pretty results:
        if (printDebugInformation) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info("\n" + prettyQueryResults);
        }

        return result;
    }

}
