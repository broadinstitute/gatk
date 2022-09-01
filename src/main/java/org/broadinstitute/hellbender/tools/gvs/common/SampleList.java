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
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryResultAndStatistics;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

public class SampleList {
    static final Logger logger = LogManager.getLogger(SampleList.class);

    private final Map<Long, String> sampleIdToNameMap = new HashMap<>();
    private long bigQueryQueryByteScanned = 0L;

    public SampleList(String sampleTableName, File sampleFile, String executionProjectId, boolean printDebugInformation, String originTool) {
        if (sampleTableName != null) {
            initializeMaps(new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS), executionProjectId, printDebugInformation, Optional.ofNullable(originTool));
        } else if (sampleFile != null) {
            initializeMaps(sampleFile);
        } else {
            throw new IllegalArgumentException("--sample-file or --sample-table must be provided.");
        }
    }

    public int size() {
        return sampleIdToNameMap.size();
    }


    public Map<Long, String> getSampleIdToNameMap() {
        return sampleIdToNameMap;
    }

    public long getBigQueryQueryByteScanned() {
        return bigQueryQueryByteScanned;
    }

    /**
     * This overload is called only from the filter creation step and should honor the `sample_info.withdrawn` column.
     */
    protected void initializeMaps(TableReference sampleTable, String executionProjectId, boolean printDebugInformation, Optional<String> originTool) {
        TableResult queryResults = querySampleTable(sampleTable.getFQTableName(), executionProjectId, printDebugInformation, originTool);

        // Add our samples to our map:
        for (final FieldValueList row : queryResults.iterateAll()) {
            long id = row.get(0).getLongValue();
            String name = row.get(1).getStringValue();
            sampleIdToNameMap.put(id, name);
        }
    }

    /**
     * This overload is called only from cohort extract and does *not* honor the `sample_info.withdrawn` column.
     */
    protected void initializeMaps(File cohortSampleFile) {
        try {
            Files.readAllLines(cohortSampleFile.toPath(), StandardCharsets.US_ASCII).stream()
                .map(s -> s.split(","))
                .forEach(tokens -> {
                    long id = Long.parseLong(tokens[0]);
                    String name = tokens[1];
                    sampleIdToNameMap.put(id, name);
                });
        } catch (IOException e) {
            throw new IllegalArgumentException("Could not parse --cohort-sample-file", e);
        }
    }

    private TableResult querySampleTable(
            String fqSampleTableName, String executionProjectId, boolean printDebugInformation, Optional<String> originTool) {
        // Get the query string:
        final String sampleListQueryString = String.format(
                "SELECT %s, %s FROM `%s` WHERE is_loaded IS TRUE and withdrawn IS NULL",
                SchemaUtils.SAMPLE_ID_FIELD_NAME,
                SchemaUtils.SAMPLE_NAME_FIELD_NAME,
                fqSampleTableName);

        Map<String, String> labelForQuery = new HashMap<>();
        if (originTool.isPresent()) {
            String originToolName = originTool.get().replaceAll("\\s","-").toLowerCase();
            labelForQuery.put("gvs_tool_name", originToolName);
            labelForQuery.put("gvs_query_name", "sample-list-creation");
        }

        // Execute the query:
        final BigQueryResultAndStatistics results = BigQueryUtils.executeQuery(
                BigQueryUtils.getBigQueryEndPoint(executionProjectId), sampleListQueryString,
                false, labelForQuery);
        bigQueryQueryByteScanned += results.queryStatistics.getTotalBytesProcessed();

        // Show our pretty results:
        if (printDebugInformation) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(results.result);
            logger.info("\n" + prettyQueryResults);
        }

        return results.result;
    }

    public static Map<Integer, LinkedList<Set<Long>>> mapSampleIdsToTableIndexes(Set<Long> sampleIds) {
        Map<Integer, LinkedList<Set<Long>>> tableMap = new HashMap<>();

        Long maxSampleId = sampleIds.stream().max(Long::compare).orElseThrow(
                () -> new GATKException("Unable to calculate max sample id, sample list may be empty")
        );

        for (Long sampleId : sampleIds) {
            int sampleTableNumber = IngestUtils.getTableNumber(sampleId, IngestConstants.partitionPerTable);

            if (!tableMap.containsKey(sampleTableNumber)) {
                tableMap.put(sampleTableNumber, new LinkedList<>());
            }

            // add the sample id to the end of the last element unless it is over
            // 1000 elements in which case make a new list
            if (tableMap.get(sampleTableNumber).size() == 0 || tableMap.get(sampleTableNumber).getLast().size() > 1000) {
                TreeSet<Long> ts = new TreeSet<>();
                tableMap.get(sampleTableNumber).addLast(ts);
            }

            tableMap.get(sampleTableNumber).getLast().add(sampleId);
        }
        return tableMap;
    }
}
