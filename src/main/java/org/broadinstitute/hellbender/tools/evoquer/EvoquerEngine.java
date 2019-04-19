package org.broadinstitute.hellbender.tools.evoquer;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.Schema;
import com.google.cloud.bigquery.TableResult;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

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
     * @param sequenceDictionary The SequenceDictionary of the reference on which the variants are based.
     * @return A {@link VCFHeader} object representing the header for all variants that have been queried from the BigQuery backend.
     */
    public VCFHeader generateVcfHeader(final Set<VCFHeaderLine> defaultHeaderLines,
                                       final SAMSequenceDictionary sequenceDictionary) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        headerLines.addAll( getEvoquerVcfHeaderLines() );
        headerLines.addAll( defaultHeaderLines );

        final VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(sequenceDictionary);

        return header;
    }

    //==================================================================================================================
    // Private Instance Methods:

    private Set<VCFHeaderLine> getEvoquerVcfHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        // Add standard VCF fields first:
        VCFStandardHeaderLines.addStandardInfoLines( headerLines, true,
                VCFConstants.STRAND_BIAS_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.STRAND_BIAS_KEY
        );

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS
        );

        // Now add GATK VCF fields:
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MAPPING_QUALITY_DEPTH));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.READ_POS_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MIN_DP_FORMAT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));

        return headerLines;
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
                // Add 1 because in the DB right now starts are exclusive:
                .start( row.get("start_position").getLongValue() +1 )
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

        final String limit_string = "LIMIT 10";

        return "SELECT " + "\n" +
                "  reference_name, start_position, end_position, reference_bases, alternate_bases, names, quality, " + "\n" +
                "  filter, call, BaseQRankSum, ClippingRankSum, variants.DP AS DP, ExcessHet, MQ, MQRankSum, MQ_DP, " + "\n" +
                "  QUALapprox, RAW_MQ, ReadPosRankSum, VarDP, variant_samples.category" + "\n" +
                "FROM " +  "\n" +
                "  `" + getFQContigTable(interval) + "` AS variant_samples " + "\n" +
                "INNER JOIN " + "\n" +
                " `" + getFQTableName(VARIANT_DATA_TABLE) + "` AS variants ON variants.end_position = variant_samples.position, " + "\n" +
                "UNNEST(variants.call) AS samples," + "\n" +
                "UNNEST(variants.alternate_bases) AS alt_bases" + "\n" +
                "WHERE " + "\n" +
                "  reference_name = '" + interval.getContig() + "' AND" + "\n" +
                "  samples.name = variant_samples.sample_id AND" + "\n" +
                "  alt_bases.alt != '<NON_REF>' AND" + "\n" +
                // Since position corresponds to end_position, we don't need to subtract 1 from thbe start here: "\n" +
                "  (position >= " + interval.getStart() + " AND position <= " + interval.getEnd() + ") AND " + "\n" +
                "  variant_samples.category = 'v' " + "\n" +
                "ORDER BY reference_name, start_position, end_position" + "\n" +
                limit_string;
    }

    //==================================================================================================================
    // Helper Data Types:

    // Example query: select count(*) from `broad-dsp-spec-ops`.gcp_joint_genotyping.variant_transforms_uuid_10 as variants, UNNEST(variants.call) as samples WHERE samples.name = "0131-01" limit 1;
    // SELECT UNNEST(call).name FROM `broad-dsp-spec-ops.gcp_joint_genotyping.variant_transforms_uuid_100`

}
