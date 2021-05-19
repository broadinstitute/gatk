package org.broadinstitute.hellbender.tools.gvs.common;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

public class SampleList {
    static final Logger logger = LogManager.getLogger(SampleList.class);

    private Map<Long, String> sampleIdMap = new HashMap<>();
    private Map<String, Long> sampleNameMap = new HashMap<>();

    public SampleList(String sampleTableName, File sampleFile, String executionProjectId, boolean printDebugInformation, String originTool) {
        if (sampleTableName != null && originTool != null) {
            initializeMaps(new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS), executionProjectId, printDebugInformation, Optional.ofNullable(originTool));
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

    protected void initializeMaps(TableReference sampleTable, String executionProjectId, boolean printDebugInformation, Optional<String> originTool) {
        TableResult queryResults = querySampleTable(sampleTable.getFQTableName(), "", executionProjectId, printDebugInformation, originTool);

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

    private TableResult querySampleTable(
        String fqSampleTableName, String whereClause, String executionProjectId, boolean printDebugInformation, Optional<String> originTool) {
        // Get the query string:
        final String sampleListQueryString =
                "SELECT " + SchemaUtils.SAMPLE_ID_FIELD_NAME + ", " + SchemaUtils.SAMPLE_NAME_FIELD_NAME +
                " FROM `" + fqSampleTableName + "`" + whereClause;

        Map<String, String> labelForQuery = new HashMap<String, String>();
        if (originTool.isPresent()) {
            String originToolName = originTool.get().toLowerCase().strip().replace(" ", "-");
            labelForQuery.put("gvs_tool_name", originToolName);
            labelForQuery.put("gvs_query_name", "sample-list-creation");
        }

        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(BigQueryUtils.getBigQueryEndPoint(executionProjectId) , sampleListQueryString, false, labelForQuery);


        // Show our pretty results:
        if (printDebugInformation) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info("\n" + prettyQueryResults);
        }

        return result;
    }

}
