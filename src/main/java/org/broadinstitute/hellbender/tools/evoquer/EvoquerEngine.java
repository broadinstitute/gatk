package org.broadinstitute.hellbender.tools.evoquer;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.Schema;
import com.google.cloud.bigquery.TableResult;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.math3.stat.descriptive.rank.Median;
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
 * Performs the work for {@link Evoquer} in a manner that is compartmentalized and reusable.
 *
 * Created by jonn on 4/17/19.
 */
class EvoquerEngine {
    private static final Logger logger = LogManager.getLogger(EvoquerEngine.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    private static final String RAW_MAPPING_QUALITY_WITH_DEPTH_KEY_SEPARATOR = ",";

    /** ID of the project containing the dataset and tables from which to pull variant data. */
    private static final String PROJECT_ID = "broad-dsp-spec-ops";

    /** ID of the table containing the names of all samples in the variant table. */
    private static final String SAMPLE_TABLE = "joint_genotyping_chr20_3_samples.sample_list";

    /**
     * The conf threshold above which variants are not included in the position tables.
     * This value is used to construct the genotype information of those missing samples
     * when they are merged together into a {@link VariantContext} object in {@link #createHighConfRefSampleGenotype(String, int, Allele, int)}.
     */
    private static final int MISSING_CONF_THRESHOLD = 60;

    /**
     * Value to insert for strand bias for the reference in high-confidence variant sample data which is missing from
     * the database.
     */
    private static final int HIGH_CONF_REFERENCE_STRAND_BIAS = 50;


    /**
     * Value to insert for MQ, ReadPosRankSum, and MQRankSum for the reference in high-confidence variant sample data which is missing from
     * the database.
     */
    private static final int MISSING_MQ_AND_READ_POS_RANK_SUM_DEFAULT_VALUE = 20;

    /**
     * Map between contig name and the BigQuery table containing position data from that contig.
     */
    private static final Map<String, String> contigPositionExpandedTableMap;

    /**
     * Map between contig name and the BigQuery table containing variant data from that contig.
     */
    private static final Map<String, String> contigVariantTableMap;

    static {
        final Map<String, String> tmpContigTableMap = new HashMap<>();
        tmpContigTableMap.put("chr20", "joint_genotyping_chr20_3_samples.pet");

        contigPositionExpandedTableMap = Collections.unmodifiableMap(tmpContigTableMap);

        final Map<String, String> tmpVariantTableMap = new HashMap<>();
        tmpVariantTableMap.put("chr20", "joint_genotyping_chr20_3_samples.vet");

        contigVariantTableMap = Collections.unmodifiableMap(tmpVariantTableMap);
    }

    //==================================================================================================================
    // Private Members:

    /** Set of sample names seen in the variant data from BigQuery. */
    private final Set<String> sampleNames = new HashSet<>();

    private final int queryRecordLimit;

    //==================================================================================================================
    // Constructors:
    EvoquerEngine(final int queryRecordLimit) {
        this.queryRecordLimit = queryRecordLimit;
    }

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
    List<VariantContext> evokeIntervals(final List<SimpleInterval> intervalList) {

        // Get the samples used in the dataset:
        populateSampleNames();

        // Now get our intervals into variants:
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
    VCFHeader generateVcfHeader(final Set<VCFHeaderLine> defaultHeaderLines,
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

        // TODO: Get a list of all possible values here so that we can make sure they're in the VCF Header!

        // Add standard VCF fields first:
        VCFStandardHeaderLines.addStandardInfoLines( headerLines, true,
                VCFConstants.STRAND_BIAS_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.RMS_MAPPING_QUALITY_KEY,
                VCFConstants.ALLELE_COUNT_KEY,
                VCFConstants.ALLELE_FREQUENCY_KEY,
                VCFConstants.ALLELE_NUMBER_KEY
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
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.FISHER_STRAND_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MIN_DP_FORMAT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));

        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        return headerLines;
    }

    /**
     * Connects to the BigQuery table for the given interval and pulls out the information on the samples that
     * contain variants.
     * @param interval {@link SimpleInterval} over which to query the BigQuery table.
     * @return A {@link List<VariantContext>} containing variants in the given {@code interval} in the BigQuery table.
     */
    private List<VariantContext> evokeInterval(final SimpleInterval interval) {

        if ( contigPositionExpandedTableMap.containsKey(interval.getContig()) ) {
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
            return createVariantsFromTableResult( result ).stream()
                    .map(GnarlyGenotyperEngine::finalizeGenotype)
                    .filter(Objects::nonNull)
                    .collect(Collectors.toList());
        }
        else {
            logger.warn("Contig missing from contigPositionExpandedTableMap, ignoring interval: " + interval.toString());
        }
        return Collections.emptyList();
    }

    private static String getPositionTableForContig(final String contig ) {
        return contigPositionExpandedTableMap.get(contig);
    }

    private static String getVariantTableForContig(final String contig ) {
        return contigVariantTableMap.get(contig);
    }

    private static String getTableQualifier() {
        return PROJECT_ID;
    }

    private static String getFQTableName( final String tableName ) {
        return getTableQualifier() + "." + tableName;
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the position
     * data specified in the given {@code interval}.
     *
     * Uses {@link #PROJECT_ID} for the project of the BigQuery table.
     * Assumes the tables have dataset information in them.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private static String getFQPositionTable(final SimpleInterval interval) {
        return getFQTableName(getPositionTableForContig( interval.getContig() ));
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the variant
     * data specified in the given {@code interval}.
     *
     * Uses {@link #PROJECT_ID} for the project of the BigQuery table.
     * Assumes the tables have dataset information in them.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private static String getFQVariantTable(final SimpleInterval interval) {
        return getFQTableName(getVariantTableForContig( interval.getContig() ));
    }

    private List<VariantContext> createVariantsFromTableResult(final TableResult result) {
        // Have to convert to int here.
        // Sloppy, but if we ever go larger than MAXINT, we have bigger problems.
        final List<VariantContext> mergedVariantContextList = new ArrayList<>((int)result.getTotalRows());

        // Details of the variants to store:
        List<VariantDetailData> currentVariantDetails = new ArrayList<>( sampleNames.size() / 10 );

        // Position / alleles are the key here.
        VariantBaseData currentVariantBaseData = new VariantBaseData();

        // Initialize with values from the first result:
        if ( result.getTotalRows() != 0 ) {
            final FieldValueList row = result.iterateAll().iterator().next();

            currentVariantBaseData.contig = row.get("reference_name").getStringValue();
            currentVariantBaseData.start = row.get("start_position").getLongValue();
            currentVariantBaseData.stop = row.get("end_position").getLongValue();

            // Fill in alleles:
            currentVariantBaseData.alleles.add( Allele.create(row.get("reference_bases").getStringValue(), true) );
            currentVariantBaseData.alleles.addAll(
                    row.get("alternate_bases").getRepeatedValue().stream()
                    .map( fieldValue -> Allele.create(fieldValue.getRecordValue().get(0).getStringValue()) )
                    .collect(Collectors.toList())
            );
        }

        for ( final FieldValueList row : result.iterateAll() ) {

            final VariantBaseData variantBaseData = new VariantBaseData();
            final VariantDetailData variantDetailData = new VariantDetailData();

            // Fill in trivial stuff:
            addBasicFieldsToVariantBuilder(row, variantBaseData);

            // Do we have a new position / new alleles?
            // If so we merge the accumulated variants and setup the new stores:
            if ( (!variantBaseData.contig.equals(currentVariantBaseData.contig)) ||
                    (variantBaseData.start != currentVariantBaseData.start)  ||
                    (variantBaseData.stop != currentVariantBaseData.stop)  ||
                    (!variantBaseData.alleles.equals(currentVariantBaseData.alleles)) ) {

                mergedVariantContextList.add(
                        mergeVariantDetails( variantBaseData, currentVariantDetails )
                );

                // Setup new values:
                currentVariantBaseData = new VariantBaseData();
                currentVariantBaseData.contig = variantBaseData.contig;
                currentVariantBaseData.start = variantBaseData.start;
                currentVariantBaseData.stop = variantBaseData.stop;
                currentVariantBaseData.alleles = variantBaseData.alleles;
                currentVariantDetails = new ArrayList<>( sampleNames.size() / 10 );
            }

            // Fill in info field stuff:
            addInfoFieldsToVariantBuilder( row, variantDetailData );

            // Fill in sample field stuff:
            // The "call" field has the genotype / sample information in it.
            // It should never be null.
            addSampleFieldsToVariantBuilder( row, result.getSchema(), variantDetailData );

            // Add our variant data to the accumulated list:
            currentVariantDetails.add( variantDetailData );
        }

        // We must merge the remaining variant details together if any are left:
        if ( !currentVariantDetails.isEmpty() ) {
            mergedVariantContextList.add(
                    mergeVariantDetails( currentVariantBaseData, currentVariantDetails )
            );
        }

        return mergedVariantContextList;
    }

    /**
     * Merges together the given variant details into a {@link VariantContext} object.
     *
     * Merges according to the following rules:
     *
     * The set of combine operations can be found here: https://github.com/Intel-HLS/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB#fields-information
     * It's also missing a new operation, which is the combined histograms (given the same double precision, combine the counts of things with the same value)
     *
     * What it lists for GATK is:
     *
     * QUAL:            set to missing
     *
     * INFO DP:         sum
     *
     * MQ:              median
     * MQRankSum:       median
     * RAW_MQ:          sum
     * QUALapprox:      sum
     * ReadPosRankSum:  median
     *
     * BaseQRankSum:    median
     * ClippingRankSum: median
     * MQ0:             median
     * ExcessHet:       median
     *
     * The GVCFs will have a FORMAT annotation for a read strand contingency table with the key "SB".
     * We combine those together across samples as element-by-element adds and move the sum to the INFO field.
     * (This gets used to calculate FS and SOR.)
     *
     * Allele-specific annotations have the same data as the traditional annotations, but the data for each
     * alternate allele is separated by a pipe as a delimiter.
     *
     * For allele-specific annotations we do a better job keeping (raw) data than the classic ones:
     * AS_RAW_*RankSum is combined by sum of histograms
     * AS_RAW_MQ is ebe sum
     * AS_SB_TABLE is ebe sum
     * AS_QD we can ignore for now.
     *
     * @param variantBaseData The {@link VariantBaseData} for the variant we are aggregating.
     * @param variantDetails {@link VariantDetailData} for each sample occuring at the given locus to be merged.
     * @return A {@link VariantContext} that combines all the information in the given {@code variantDetails}.
     */
    private VariantContext mergeVariantDetails(final VariantBaseData variantBaseData,
                                               final List<VariantDetailData> variantDetails) {

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

        // Populate trivial fields in variant context builder:
        variantContextBuilder.chr(variantBaseData.contig)
                            .start(variantBaseData.start)
                            .stop(variantBaseData.stop)
                            .alleles(variantBaseData.alleles);

        // Get the alleles we actually care about here:
        final List<Allele> genotypeAlleles = variantBaseData.alleles.stream().filter( a -> !a.equals(Allele.NON_REF_ALLELE) ).collect(Collectors.toList());

        // no need to populate ID
        // no need to populate FILTER
        // no need to populate QUAL as per rules

        final Set<String>    samplesMissing  = new HashSet<>( sampleNames );
        final List<Genotype> sampleGenotypes = new ArrayList<>( sampleNames.size() );

        int depth = 0;
        int variantDepth = 0;
        int mapQualityDepth = 0;
        long rawMq = 0;
        double qualApprox = 0;

        final List<Double> mq             = new ArrayList<>(variantDetails.size());
        final List<Double> mqRankSum      = new ArrayList<>(variantDetails.size());
        final List<Double> readPosRankSum = new ArrayList<>(variantDetails.size());

        // Now go over each sample and aggregate the data as per the rules:
        for ( final VariantDetailData sampleData : variantDetails  ) {

            // ------------------------------------------------------------------
            // INFO fields:

            // Simple aggregations on:
            //   DP, MQ_DP, QUALapprox, RAW_MQ, VAR_DP
            if ( sampleData.infoDp != null ) {
                depth += sampleData.infoDp;
            }
            if ( sampleData.infoQualApprox != null ) {
                qualApprox += sampleData.infoQualApprox;
            }
            if ( sampleData.infoRawMq != null ) {
                rawMq += sampleData.infoRawMq;
            }
            if ( sampleData.infoMqDp != null ) {
                // TODO: This may not be right!
                mapQualityDepth += sampleData.infoMqDp;
            }
            if ( sampleData.infoVarDp != null ) {
                // TODO: This may not be right!
                variantDepth += sampleData.infoVarDp;
            }

            // Median calculations on:
            //   MQ, MQRankSum, readPosRankSum

            if ( sampleData.infoMq != null ) {
                mq.add(sampleData.infoMq);
            }
            if ( sampleData.infoMqRankSum != null ) {
                mqRankSum.add(sampleData.infoMqRankSum);
            }
            if ( sampleData.infoReadPosRankSum != null ) {
                readPosRankSum.add(sampleData.infoReadPosRankSum);
            }

            // Genotype fields should just be added as-is:
            sampleGenotypes.add( createGenotypeFromVariantDetails(sampleData, genotypeAlleles) );

            // Account for this sample in our sample set:
            samplesMissing.remove( sampleData.sampleName );
        }

        // Calculate the median values:
        double mqMedian             = new Median().evaluate(mq.stream().mapToDouble(Double::doubleValue).toArray());
        double mqRankSumMedian      = new Median().evaluate(mqRankSum.stream().mapToDouble(Double::doubleValue).toArray());
        double readPosRankSumMedian = new Median().evaluate(readPosRankSum.stream().mapToDouble(Double::doubleValue).toArray());

        // Ensure no NaN values:
        // TODO: Make sure these values are good:
        mqMedian             = Double.isNaN(mqMedian) ? MISSING_MQ_AND_READ_POS_RANK_SUM_DEFAULT_VALUE : mqMedian;
        mqRankSumMedian      = Double.isNaN(mqRankSumMedian) ? MISSING_MQ_AND_READ_POS_RANK_SUM_DEFAULT_VALUE : mqRankSumMedian;
        readPosRankSumMedian = Double.isNaN(readPosRankSumMedian) ? MISSING_MQ_AND_READ_POS_RANK_SUM_DEFAULT_VALUE : readPosRankSumMedian;

        // Add the aggregate values to our builder:
        variantContextBuilder.attribute(VCFConstants.DEPTH_KEY, depth);
        variantContextBuilder.attribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY, qualApprox);
        variantContextBuilder.attribute(GATKVCFConstants.MAPPING_QUALITY_DEPTH, mapQualityDepth);
        variantContextBuilder.attribute(GATKVCFConstants.VARIANT_DEPTH_KEY, variantDepth);
        variantContextBuilder.attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, String.format("%d%s%d", rawMq, RAW_MAPPING_QUALITY_WITH_DEPTH_KEY_SEPARATOR, depth));

        variantContextBuilder.attribute(VCFConstants.RMS_MAPPING_QUALITY_KEY, mqMedian);
        variantContextBuilder.attribute(GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY, mqRankSumMedian);
        variantContextBuilder.attribute(GATKVCFConstants.READ_POS_RANK_SUM_KEY, readPosRankSumMedian);

        // Now add in empty values for each sample that was not in our variant details.
        // We assume these samples have high confidence reference regions at this allele
        for ( final String sample : samplesMissing ) {
            sampleGenotypes.add( createHighConfRefSampleGenotype(sample, depth, genotypeAlleles.get(0), genotypeAlleles.size() ) );
        }

        // Set our genotypes:
        variantContextBuilder.genotypes( sampleGenotypes );

        // Return the VC:
        return variantContextBuilder.make();
    }

    /**
     * Create a {@link Genotype} from the given {@link VariantDetailData} and {@code alleles.}
     * @param sampleData {@link VariantDetailData} containing sample information from which to create a {@link Genotype}.
     * @param alleles {@link Allele} objects from which to create the {@link Genotype}.
     * @return A {@link Genotype} object containing the information in the given {@code sampleData}.
     */
    private Genotype createGenotypeFromVariantDetails(final VariantDetailData sampleData, final List<Allele> alleles) {
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
        genotypeBuilder.name(sampleData.sampleName);

        // GT:AD:DP:GQ:PL:SB
        genotypeBuilder.alleles( alleles );
        genotypeBuilder.AD( sampleData.gtAd );
        genotypeBuilder.DP( sampleData.gtDp );
        genotypeBuilder.GQ( sampleData.gtGq );
        genotypeBuilder.PL( sampleData.gtPl );
        genotypeBuilder.attribute(
                VCFConstants.STRAND_BIAS_KEY,
                Arrays.stream(sampleData.gtSb).map(i -> Integer.toString(i)).collect(Collectors.joining(","))
        );
        return genotypeBuilder.make();
    }

    /**
     * Creates a {@link Genotype} object containing default "high-confidence" field values for the given
     * {@code sampleId}.
     *
     * These values are placeholders used only so the resuling {@link VariantContext} will contain information from the
     * given {@code sampleId} to enable genotyping.  These data are not necessarily reflective of the actual sample
     * information.
     *
     * @param sampleId The ID of a sample for which to generate default information in the given {@link VariantContextBuilder}.
     * @param depth The depth to use for the sample.
     * @param refAllele The reference {@link Allele} in this variant.
     * @param numAlleles Total number of alleles in this variant, including the reference allele.
     * @return The {@link Genotype} object with default "high-confidence" field values for the given {@code sampleId}.
     */
    private Genotype createHighConfRefSampleGenotype(final String sampleId,
                                                     final int depth,
                                                     final Allele refAllele,
                                                     final int numAlleles ) {

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        genotypeBuilder.name(sampleId);

        // GT:AD:DP:GQ:PL:SB
        // TODO: this is probably wrong - need to remove one copy for the <NON_REF> allele so we can ignore it.
        genotypeBuilder.alleles( Collections.nCopies( numAlleles, refAllele ) );
        genotypeBuilder.AD( Collections.nCopies( numAlleles, depth ).stream().mapToInt( i -> i).toArray() );
        genotypeBuilder.DP(depth);
        genotypeBuilder.GQ(MISSING_CONF_THRESHOLD);

        // Setup our PLs:
        final List<Integer> pls = new ArrayList<>((numAlleles-1)*3);
        pls.add(0);
        for ( int i = 0 ; i < ((numAlleles-1)*3)-1 ; i++) { pls.add(MISSING_CONF_THRESHOLD); }
        genotypeBuilder.PL( pls.stream().mapToInt(i -> i).toArray() );

        // Setup our SBs:
        final List<Integer> sbs = new ArrayList<>(numAlleles * 2);
        sbs.add(HIGH_CONF_REFERENCE_STRAND_BIAS);
        sbs.add(HIGH_CONF_REFERENCE_STRAND_BIAS);
        for ( int i = 0 ; i < (numAlleles-1)*2 ; i++) { sbs.add(0); }
        genotypeBuilder.attribute(
                VCFConstants.STRAND_BIAS_KEY,
                sbs.stream().map( Object::toString ).collect(Collectors.joining(","))
        );

        return genotypeBuilder.make();
    }

    private void addBasicFieldsToVariantBuilder(final FieldValueList row,
                                                final VariantBaseData variantBaseData) {

        variantBaseData.contig = row.get("reference_name").getStringValue();
        variantBaseData.start = row.get("start_position").getLongValue();
        variantBaseData.stop = row.get("end_position").getLongValue();

        // Get the filter(s):
        if ( !row.get("filter").isNull() ) {
            variantBaseData.filters =
                    row.get("filter").getRepeatedValue().stream()
                        .map( FieldValue::getStringValue )
                        .collect(Collectors.toSet());
        }

        // Fill in alleles:
        final List<Allele> alleles = new ArrayList<>(5);
        alleles.add( Allele.create(row.get("reference_bases").getStringValue(), true) );

        alleles.addAll(
                row.get("alternate_bases").getRepeatedValue().stream()
                    .map( fieldValue -> Allele.create(fieldValue.getRecordValue().get(0).getStringValue()) )
                    .collect(Collectors.toList())
        );

        // Add the alleles:
        variantBaseData.alleles = alleles;
    }

    private void addInfoFieldsToVariantBuilder(final FieldValueList row,
                                               final VariantDetailData variantDetailData) {

        // NOTE: If a field is null we ignore it for the aggregation step (as per Laura's instructions).
        if ( !row.get(VCFConstants.DEPTH_KEY).isNull() ) {
            variantDetailData.infoDp = (int) row.get(VCFConstants.DEPTH_KEY).getLongValue();
        }
        if ( !row.get(VCFConstants.RMS_MAPPING_QUALITY_KEY).isNull() ) {
            variantDetailData.infoMq = row.get(VCFConstants.RMS_MAPPING_QUALITY_KEY).getDoubleValue();
        }
        if ( !row.get(GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY).isNull() ) {
            variantDetailData.infoMqRankSum = row.get(GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY).getDoubleValue();
        }
        if ( !row.get(GATKVCFConstants.MAPPING_QUALITY_DEPTH).isNull() ) {
            variantDetailData.infoMqDp = (int) row.get(GATKVCFConstants.MAPPING_QUALITY_DEPTH).getLongValue();
        }
        if ( !row.get(GATKVCFConstants.RAW_QUAL_APPROX_KEY).isNull() ) {
            variantDetailData.infoQualApprox = (int) row.get(GATKVCFConstants.RAW_QUAL_APPROX_KEY).getLongValue();
        }
        if ( !row.get(GATKVCFConstants.READ_POS_RANK_SUM_KEY).isNull() ) {
            variantDetailData.infoReadPosRankSum = row.get(GATKVCFConstants.READ_POS_RANK_SUM_KEY).getDoubleValue();
        }
        if ( !row.get(GATKVCFConstants.VARIANT_DEPTH_KEY).isNull() ) {
            variantDetailData.infoVarDp = (int) row.get(GATKVCFConstants.VARIANT_DEPTH_KEY).getLongValue();
        }
        if ( !row.get(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY).isNull() ) {
            variantDetailData.infoRawMq = Math.round(row.get(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY).getDoubleValue());
        }
    }

    private void addSampleFieldsToVariantBuilder( final FieldValueList row,
                                                  final Schema schema,
                                                  final VariantDetailData variantDetailData) {

        // Create our call data with the schema information
        // to enable name-based access to fields:
        final FieldValueList callData =
                FieldValueList.of(
                        // I'm not sure why we need this extra layer of unwrapping:
                        row.get("call").getRecordValue().get(0).getRecordValue(),
                        schema.getFields().get("call").getSubFields()
                );

        final String sampleName = callData.get("name").getStringValue();

        final int dp = (int)callData.get(VCFConstants.DEPTH_KEY).getLongValue();
        final int gq = (int)callData.get(VCFConstants.GENOTYPE_QUALITY_KEY).getLongValue();

        variantDetailData.sampleName = sampleName;
        variantDetailData.gtDp = dp;
        variantDetailData.gtGq = gq;
        //TODO: Do we need to track phaseset here?
        //TODO: Do we need to track MIN_DP here?
        //TODO: Do we need to track PGT here?
        //TODO: Do we need to track PID here?

        // Get the array fields:

        // Add the alleles to our variantDetailData:
        variantDetailData.gtGenotype = callData.get("genotype").getRepeatedValue().stream()
                .mapToInt( f -> (int)f.getLongValue() )
                .toArray();

        // AD should never be null:
        variantDetailData.gtAd = callData.get(VCFConstants.GENOTYPE_ALLELE_DEPTHS).getRecordValue().stream()
                .map( FieldValue::getLongValue )
                .mapToInt( Long::intValue )
                .toArray();

        // PL should never be null:
        variantDetailData.gtPl = callData.get(VCFConstants.GENOTYPE_PL_KEY).getRecordValue().stream()
                                 .map( FieldValue::getLongValue )
                                 .mapToInt( Long::intValue )
                                 .toArray();

        if ( !callData.get(VCFConstants.STRAND_BIAS_KEY).isNull() ) {

            variantDetailData.gtSb = callData.get(VCFConstants.STRAND_BIAS_KEY).getRecordValue().stream()
                    .map( FieldValue::getLongValue )
                    .mapToInt( Long::intValue )
                    .boxed()
                    .toArray(Integer[]::new);
        }
    }

    private void populateSampleNames() {
        // Get the query string:
        final String sampleListQueryString = getSampleListQueryString();

        logger.info("Created Query: \n" + sampleListQueryString);

        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        // Show our pretty results:
        logger.info("Pretty Query Results:");
        final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
        logger.info( "\n" + prettyQueryResults );

        // Add our samples to our map:
        for ( final FieldValueList row : result.iterateAll() ) {
            sampleNames.add( row.get(0).getStringValue() );
        }
    }

    private String getVariantQueryString( final SimpleInterval interval ) {

        String limitString = "";
        if ( queryRecordLimit > 0 ) {
            limitString = "LIMIT " + queryRecordLimit;
        }

        // TODO: Replace column names with variable field values for consistency (as above)
        return "SELECT " + "\n" +
                "  reference_name, start_position, end_position, reference_bases, alternate_bases, names, quality, " + "\n" +
                "  filter, call, BaseQRankSum, ClippingRankSum, variants.DP AS DP, ExcessHet, MQ, " + "\n" +
                "  MQRankSum, MQ_DP, QUALapprox, RAW_MQ, ReadPosRankSum, VarDP, variant_samples.state" + "\n" +
                "FROM " +  "\n" +
                "  `" + getFQPositionTable(interval) + "` AS variant_samples " + "\n" +
                "INNER JOIN " + "\n" +
                " `" + getFQVariantTable(interval) + "` AS variants ON variants.end_position = variant_samples.position, " + "\n" +
                "UNNEST(variants.call) AS samples" + "\n" +
                "WHERE " + "\n" +
                "  reference_name = '" + interval.getContig() + "' AND" + "\n" +
                "  samples.name = variant_samples.sample AND" + "\n" +
                // Since position corresponds to end_position, we don't need to subtract 1 from thbe start here: "\n" +
                "  (position >= " + interval.getStart() + " AND position <= " + interval.getEnd() + ") AND " + "\n" +
                "  variant_samples.state = 'v' " + "\n" +
                "ORDER BY reference_name, start_position, end_position" + "\n" +
                limitString;
    }

    private static String getSampleListQueryString() {
        return "SELECT sample FROM `" + getFQTableName(SAMPLE_TABLE)+ "`";
    }

    //==================================================================================================================
    // Helper Data Types:

    /**
     * A class to hold high-level variant information without constructing an entire {@link VariantContext} object.
     */
    private class VariantBaseData {
        String contig;
        long start;
        long stop;

        List<Allele> alleles = new ArrayList<>();
        Set<String> filters = new HashSet<>();
    }

    /**
     * A class to hold variant detail information about a particular sample
     * without constructing an entire {@link VariantContext} object.
     */
    private class VariantDetailData {

        // High-level fields:
        String sampleName;

        // Info Fields:
        //DP=7;MQ=34.15;MQRankSum=1.300;MQ_DP=7;QUALapprox=10;RAW_MQandDP=239.05,7;ReadPosRankSum=1.754;VarDP=7
        Integer infoDp;
        Double infoMq;
        Double infoMqRankSum;
        Integer infoMqDp;
        Integer infoQualApprox;
        Long infoRawMq;
        Double infoReadPosRankSum;
        Integer infoVarDp;

        // Genotype fields:
        // GT:AD:DP:GQ:PL:SB 0/1:   4,3,0:  7: 10:  10,0,119,101,128,229:0,4,0,3
        int gtGenotype[];
        int gtAd[];
        Integer gtDp;
        Integer gtGq;
        int gtPl[];
        Integer gtSb[];
    }
}
