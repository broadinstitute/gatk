package org.broadinstitute.hellbender.tools.variantdb.arrays.tables;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SampleList {
    static final Logger logger = LogManager.getLogger(SampleList.class);

    public static final List<String> SAMPLE_LIST_FIELDS = Arrays.asList(SchemaUtils.SAMPLE_NAME_FIELD_NAME, SchemaUtils.SAMPLE_ID_FIELD_NAME);

    public static Map<String, Integer> getSampleNameMap(TableReference sampleTable, List<String> samples, boolean printDebugInformation) {
        Map<String, Integer> results = new HashMap<>();
        // create optional where clause
        String whereClause = "";
        if (samples != null && samples.size() > 0) {
            whereClause = " WHERE " + SchemaUtils.SAMPLE_NAME_FIELD_NAME + " in (\'" + StringUtils.join(samples, "\',\'") + "\') ";
        }

        TableResult queryResults = querySampleTable(sampleTable.getFQTableName(), whereClause, printDebugInformation);

        // Add our samples to our map:
        for (final FieldValueList row : queryResults.iterateAll()) {
            results.put(row.get(1).getStringValue(), (int) row.get(0).getLongValue());
        }
        return results;
    }


    public static Map<Integer, String> getSampleIdMap(TableReference sampleTable, boolean printDebugInformation) {

        Map<Integer, String> results = new HashMap<>();
        TableResult queryResults = querySampleTable(sampleTable.getFQTableName(), "", printDebugInformation);

        // Add our samples to our map:
        for (final FieldValueList row : queryResults.iterateAll()) {
            results.put((int) row.get(0).getLongValue(), row.get(1).getStringValue());
        }
        return results;
    }

    private static TableResult querySampleTable(String fqSampleTableName, String whereClause, boolean printDebugInformation) {
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
