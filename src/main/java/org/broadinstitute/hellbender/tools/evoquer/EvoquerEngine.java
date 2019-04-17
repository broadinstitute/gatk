package org.broadinstitute.hellbender.tools.evoquer;

import com.google.cloud.bigquery.TableResult;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * EvoquerEngine ("EvokerEngine"):
 * Extract Variants Out of big QUERy Engine.
 *
 * Performs the work for {@link Evoquer} in a manner that is extensible.
 *
 * Created by jonn on 4/17/19.
 */
public class EvoquerEngine {
    private static final Logger logger = LogManager.getLogger(EvoquerEngine.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    /** ID of the project containing the dataset and tables from which to pull variant data. */
    private static final String PROJECT_ID = "broad-dsp-spec-ops";

    /** ID dataset and tables from which to pull variant data. */
    private static final String DATASET_ID = "gcp_joint_genotyping";

    /** ID dataset and tables from which to pull variant data. */
    private static final String VARIANT_DATA_TABLE = "variant_transforms_uuid_10";

    /**
     * Map between contig name and the BigQuery table containing data from that contig.
     */
    private static final Map<String, String> contigTableMap;

    static {
        final Map<String, String> tmpContigTableMap = new HashMap<>();
        tmpContigTableMap.put("chr2", "chr2_sample_100_new_way");

        contigTableMap = Collections.unmodifiableMap(tmpContigTableMap);
    }

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    public EvoquerEngine() {}

    //==================================================================================================================
    // Override Methods:

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Public Instance Methods:

    /**
     * Connects to the BigQuery table for the given {@link List<SimpleInterval>} and pulls out the information on the samples that
     * contain variants.
     * @param intervalList {@link List<SimpleInterval>} over which to query the BigQuery table.
     * @return A {@link List<VariantContext>} containing variants in the given {@code interval} in the BigQuery table.
     */
    public List<VariantContext> evokeIntervals(final List<SimpleInterval> intervalList) {

        return intervalList.stream()
                .flatMap( interval -> evokeInterval(interval).stream() )
                .collect(Collectors.toList());
    }

    //==================================================================================================================
    // Private Instance Methods:

    /**
     * Connects to the BigQuery table for the given interval and pulls out the information on the samples that
     * contain variants.
     * @param interval {@link SimpleInterval} over which to query the BigQuery table.
     * @return A {@link List<VariantContext>} containing variants in the given {@code interval} in the BigQuery table.
     */
    private List<VariantContext> evokeInterval(final SimpleInterval interval) {

        if ( contigTableMap.containsKey(interval.getContig()) ) {
            // Get the query string:
            final String variantQueryString = getVariantQueryString(interval);

            logger.info("Created Query: \n" + variantQueryString);

            // Execute the query:
            final TableResult result = BigQueryUtils.executeQuery(variantQueryString);

            // Show our pretty results:
            logger.info("Pretty Query Results:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info( "\n" + prettyQueryResults );
        }
        else {
            logger.warn("Contig missing from contigTableMap, ignoring interval: " + interval.toString());
        }
        return Collections.emptyList();
    }

    private static String getTableForContig( final String contig ) {
        return contigTableMap.get(contig);
    }

    private static String getTableQualifier() {
        return PROJECT_ID + "." + DATASET_ID;
    }

    private static String getFQTableName( final String tableName ) {
        return getTableQualifier() + "." + tableName;
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the data specified
     * in the given {@code interval}.
     *
     * Uses {@link #PROJECT_ID} and {@link #DATASET_ID} for the project and dataset of the BigQuery table.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private static String getFQContigTable(final SimpleInterval interval) {
        return getFQTableName(getTableForContig( interval.getContig() ));
    }

    private String getVariantQueryString( final SimpleInterval interval ) {

        final String lim1 = "LIMIT 5";
        final String lim2 = "LIMIT 10";

        return "WITH variant_samples AS (" + "\n" +
                "  SELECT sample_id, position FROM `" + getFQContigTable(interval) + "` " + "\n" +
                "  WHERE " + "\n" +
                "    (position >= " + interval.getStart() + " AND position <= " + interval.getEnd() + ") AND " + "\n" +
                "    category = 'v' " + "\n" +
                "  " + lim1 + "\n" +
                ")" + "\n" +
                "" + "\n" +
                "SELECT " + "\n" +
                "  reference_name, start_position, end_position, reference_bases, alternate_bases, names, quality, " + "\n" +
                "  filter, call, BaseQRankSum, ClippingRankSum, variants.DP, ExcessHet, MQ, MQRankSum, MQ_DP, " + "\n" +
                "  QUALapprox, RAW_MQ, ReadPosRankSum, VarDP" + "\n" +
                "FROM" + "\n" +
                "  `" + getFQTableName(VARIANT_DATA_TABLE) + "` as variants" + "\n" +
                "CROSS JOIN UNNEST(variants.call) AS sample " + "\n" +
                "CROSS JOIN UNNEST(alternate_bases) AS alt_allele " + "\n" +
                "WHERE" + "\n" +
                "  reference_name = '" + interval.getContig() + "' AND" + "\n" +
                "  alt_allele.alt != '<NON_REF>' AND " + "\n" +
                "  sample.name IN (SELECT sample_id FROM variant_samples) AND" + "\n" +
                "  (start_position IN (SELECT position FROM variant_samples) OR" + "\n" +
                "   end_position   IN (SELECT position FROM variant_samples)) " + "\n" +
                "  " + lim2;
    }

    //==================================================================================================================
    // Helper Data Types:

    // Example query: select count(*) from `broad-dsp-spec-ops`.gcp_joint_genotyping.variant_transforms_uuid_10 as variants, UNNEST(variants.call) as samples samples.name = "0131-01" limit 1;
    // SELECT UNNEST(call).name FROM `broad-dsp-spec-ops.gcp_joint_genotyping.variant_transforms_uuid_100`

}
