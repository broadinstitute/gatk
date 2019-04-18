package org.broadinstitute.hellbender.tools.evoquer;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.Schema;
import com.google.cloud.bigquery.TableResult;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;

import java.util.*;
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

    /** Set of sample names seen in the variant data from BigQuery. */
    private Set<String> sampleNames = new HashSet<>();

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

    /**
     * Generates a {@link VCFHeader} object based on the VariantContext objects queried from the BigQuery backend.
     * If no objects have been queried, this will return a default {@link VCFHeader}.
     * @param defaultHeaderLines The default header lines to be added to the top of the VCF header.
     * @return A {@link VCFHeader} object representing the header for all variants that have been queried from the BigQuery backend.
     */
    public VCFHeader generateVcfHeader(final Set<VCFHeaderLine> defaultHeaderLines ) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        headerLines.addAll( getEvoquerVcfHeaderLines() );
        headerLines.addAll( defaultHeaderLines );

        logger.info( "VCF Header Lines:" );
        for ( final VCFHeaderLine line : headerLines ) {
            logger.info(line);
        }

        return new VCFHeader(headerLines, sampleNames);
    }

    //==================================================================================================================
    // Private Instance Methods:

    private List<VCFHeaderLine> getEvoquerVcfHeaderLines() {
        // TODO: FINISHME!
        return Collections.emptyList();
    }

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

            // Convert results into variant context objects:
            return createVariantsFromTableResult( result );
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

    private List<VariantContext> createVariantsFromTableResult(final TableResult result) {
        // Have to convert to int here.
        // Sloppy, but if we ever go larger than MAXINT, we have bigger problems.
        final List<VariantContext> variantContextList = new ArrayList<>((int)result.getTotalRows());

        for ( final FieldValueList row : result.iterateAll() ) {
            final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

            // TODO: add the info fields / genotypes / sample information!

            // Fill in trivial stuff:
            addBasicFieldsToVariantBuilder(row, variantContextBuilder);

            // Fill in info field stuff:
            addInfoFieldsToVariantBuilder( row, variantContextBuilder );

            // Fill in sample field stuff:
            // The "call" field has the genotype / sample information in it.
            // It should never be null.
            addSampleFieldsToVariantBuilder( row, result.getSchema(), variantContextBuilder );

            // Add our variant context to the list:
            variantContextList.add( variantContextBuilder.make() );
        }

        return variantContextList;
    }

    private void addBasicFieldsToVariantBuilder(final FieldValueList row,
                                                final VariantContextBuilder variantContextBuilder) {
        variantContextBuilder
                .chr( row.get("reference_name").getStringValue() )
                .start( row.get("start_position").getLongValue() )
                .stop( row.get("end_position").getLongValue() );

        // Get the filter(s):
        if ( !row.get("filter").isNull() ) {
            variantContextBuilder.filters(
                    row.get("filter").getRepeatedValue().stream()
                        .map( fieldValue -> fieldValue.getRecordValue().get(0).getStringValue() )
                        .collect(Collectors.toSet())
            );
        }

        // Qual:
        if ( !row.get("quality").isNull() ) {
            variantContextBuilder.log10PError( row.get("quality").getDoubleValue() / -10.0 );
        }

        // Fill in alleles:
        final List<String> alleles = new ArrayList<>(5);
        alleles.add( row.get("reference_bases").getStringValue() );

        alleles.addAll(
                row.get("alternate_bases").getRepeatedValue().stream()
                    .map( fieldValue -> fieldValue.getRecordValue().get(0).getStringValue() )
                    .collect(Collectors.toList())
        );

        // Add the alleles:
        variantContextBuilder.alleles( alleles );
    }

    private void addInfoFieldsToVariantBuilder(final FieldValueList row,
                                               final VariantContextBuilder variantContextBuilder) {

        // Start by adding in the END to the INFO field:
        addInfoFieldToVariantContextBuilder( row, "DP", "DP", variantContextBuilder );
        addInfoFieldToVariantContextBuilder( row, "MQ", "MQ", variantContextBuilder );
        addInfoFieldToVariantContextBuilder( row, "MQRankSum", "MQRankSum", variantContextBuilder );
        addInfoFieldToVariantContextBuilder( row, "MQ_DP", "MQ_DP", variantContextBuilder );
        addInfoFieldToVariantContextBuilder( row, "QUALapprox", "QUALapprox", variantContextBuilder );
        addInfoFieldToVariantContextBuilder( row, "ReadPosRankSum", "ReadPosRankSum", variantContextBuilder );
        addInfoFieldToVariantContextBuilder( row, "VarDP", "VarDP", variantContextBuilder );

        // Handle RAW_MQandDP field:
        final String dp = row.get("DP").getStringValue();
        final String raw_mq = row.get("RAW_MQ").getStringValue();
        variantContextBuilder.attribute("RAW_MQandDP", raw_mq + "," + dp);
    }

    private void addSampleFieldsToVariantBuilder( final FieldValueList row,
                                                  final Schema schema,
                                                  final VariantContextBuilder variantContextBuilder) {

        // Create our call data with the schema information
        // to enable name-based access to fields:
        final FieldValueList callData =
                FieldValueList.of(
                        // I'm not sure why we need this extra layer of unwrapping:
                        row.get("call").getRecordValue().get(0).getRecordValue(),
                        schema.getFields().get("call").getSubFields()
                );

        final Allele refAllele = Allele.create(row.get("reference_bases").getStringValue().getBytes(), true);

        final String sampleName = callData.get("name").getStringValue();
        final int dp = (int)callData.get("DP").getLongValue();
        final int gq = (int)callData.get("GQ").getLongValue();

        // Handle the ref allele genotype first:
        final GenotypeBuilder refGenotypeBuilder = new GenotypeBuilder();
        refGenotypeBuilder.name(sampleName);
        refGenotypeBuilder.DP(dp);
        refGenotypeBuilder.GQ(gq);

        // Now handle the alt allele genotypes:

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        // Get the scalar fields first:
        genotypeBuilder.name(sampleName);
        genotypeBuilder.DP(dp);
        genotypeBuilder.GQ(gq);
        addScalarAttributeToGenotypeBuilder( callData, "phaseset", genotypeBuilder );
        addScalarAttributeToGenotypeBuilder( callData, "MIN_DP", genotypeBuilder );
        addScalarAttributeToGenotypeBuilder( callData, "PGT", genotypeBuilder );
        addScalarAttributeToGenotypeBuilder( callData, "PID", genotypeBuilder );

        // Get the array fields:

        // Add the alleles:
        final List<Allele> alleleList = new ArrayList<>();
        callData.get("genotype").getRepeatedValue().stream()
                .map( f -> (int)f.getLongValue() )
                .map( gtIndex ->
                    {
                        if ( gtIndex == 0 ) {
                            alleleList.add(refAllele);
                        }
                        else {
                            alleleList.add(
                                    Allele.create(
                                            row.get("alternate_bases")
                                                    .getRecordValue()
                                                    .get("alt")
                                                    .getRepeatedValue()
                                                    .get(gtIndex)
                                                    .getStringValue()
                                    )
                            );
                        }
                        return true;
                    }
                );
        genotypeBuilder.alleles( alleleList );

        // AD should not be null:
        genotypeBuilder.AD(
                callData.get("AD").getRecordValue().stream()
                        .map( FieldValue::getLongValue )
                        .mapToInt( Long::intValue )
                        .toArray()
        );

        // PL should not be null:
        genotypeBuilder.PL(
                callData.get("PL").getRecordValue().stream()
                        .map( FieldValue::getLongValue )
                        .mapToInt( Long::intValue )
                        .toArray()
        );

        if ( !callData.get("SB").isNull() ) {
            genotypeBuilder.attribute("SB",
                    callData.get("SB").getRecordValue().stream()
                            .map(FieldValue::getStringValue)
                            .collect(Collectors.joining(","))
            );
        }

        variantContextBuilder.genotypes( genotypeBuilder.make() );

        // Before we leave, let's add the sample name to sampleNames:
        sampleNames.add( sampleName );
    }

    private void addScalarAttributeToGenotypeBuilder(final FieldValueList row,
                                                     final String fieldName,
                                                     final GenotypeBuilder genotypeBuilder ) {
        // Only add the info if it is not null:
        if ( !row.get(fieldName).isNull() ) {
            genotypeBuilder.attribute(fieldName, row.get(fieldName).getStringValue());
        }
    }

    private void addInfoFieldToVariantContextBuilder(final FieldValueList row,
                                                     final String columnName,
                                                     final String infoFieldName,
                                                     final VariantContextBuilder variantContextBuilder ) {
        // Only add the info if it is not null:
        if ( !row.get(columnName).isNull() ) {
            variantContextBuilder.attribute(infoFieldName, row.get(columnName).getStringValue());
        }
    }

    private String getVariantQueryString( final SimpleInterval interval ) {

        final String limit = "LIMIT 10";

        return "WITH variant_samples AS (" + "\n" +
                "  SELECT sample_id, position FROM `" + getFQContigTable(interval) + "` " + "\n" +
                "  WHERE " + "\n" +
                "    (position >= " + interval.getStart() + " AND position <= " + interval.getEnd() + ") AND " + "\n" +
                "    category = 'v' " + "\n" +
                ")" + "\n" +
                "SELECT " + "\n" +
                "  reference_name, start_position, end_position, reference_bases, alternate_bases, names, quality," +
                "  filter, call, BaseQRankSum, ClippingRankSum, variants.DP AS DP, ExcessHet, MQ, MQRankSum, MQ_DP," +
                "  QUALapprox, RAW_MQ, ReadPosRankSum, VarDP " + "\n" +
                "FROM " + "\n" +
                "  `" + getFQTableName(VARIANT_DATA_TABLE) + "` AS variants, " + "\n" +
                "UNNEST(variants.call) AS samples," + "\n" +
                "UNNEST(variants.alternate_bases) as alt_bases" + "\n" +
                "INNER JOIN variant_samples on (variants.end_position = variant_samples.position OR variants.start_position = variant_samples.position)" + "\n" +
                "WHERE " + "\n" +
                "  reference_name = '" + interval.getContig() + "' AND" + "\n" +
                "  samples.name = variant_samples.sample_id AND" + "\n" +
                "  alt_bases.alt != '<NON_REF>'" + "\n" +
                "ORDER BY reference_name, start_position, end_position" + "\n" +
                limit;
    }

    //==================================================================================================================
    // Helper Data Types:

    // Example query: select count(*) from `broad-dsp-spec-ops`.gcp_joint_genotyping.variant_transforms_uuid_10 as variants, UNNEST(variants.call) as samples WHERE samples.name = "0131-01" limit 1;
    // SELECT UNNEST(call).name FROM `broad-dsp-spec-ops.gcp_joint_genotyping.variant_transforms_uuid_100`

}
