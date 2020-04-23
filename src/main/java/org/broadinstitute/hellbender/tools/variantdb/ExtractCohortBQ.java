package org.broadinstitute.hellbender.tools.variantdb;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.util.HashSet;
import java.util.Set;

public class ExtractCohortBQ {
    private static final Logger logger = LogManager.getLogger(ExtractCohortBQ.class);

    public static Set<String> populateSampleNames(TableReference sampleTableRef, boolean printDebugInformation) {
        Set<String> results = new HashSet<>();

        // Get the query string:
        final String sampleListQueryString = "SELECT sample FROM `" + sampleTableRef.getFQTableName() + "`";
        ;

        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        // Show our pretty results:
        if (printDebugInformation) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info("\n" + prettyQueryResults);
        }

        // Add our samples to our map:
        for (final FieldValueList row : result.iterateAll()) {
            results.add(row.get(0).getStringValue());
        }

        return results;
    }

}
