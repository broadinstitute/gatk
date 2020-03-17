package org.broadinstitute.hellbender.tools.evoquer;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import com.google.errorprone.annotations.Var;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyperEngine;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

import static htsjdk.variant.vcf.VCFConstants.PASSES_FILTERS_v4;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.*;

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

    public static final String POSITION_FIELD_NAME = "position";
    public static final String VALUES_ARRAY_FIELD_NAME = "values";

    public static final String SAMPLE_FIELD_NAME = "sample";
    public static final String STATE_FIELD_NAME = "state";
    public static final String REF_ALLELE_FIELD_NAME = "ref";
    public static final String ALT_ALLELE_FIELD_NAME = "alt";

    public static final String GENOTYPE_FIELD_PREFIX = "call_";
    public static final String MULTIVALUE_FIELD_DELIMITER = ",";

    public static final ImmutableSet<String> REQUIRED_FIELDS = ImmutableSet.of(
            SAMPLE_FIELD_NAME,
            STATE_FIELD_NAME,
            REF_ALLELE_FIELD_NAME,
            ALT_ALLELE_FIELD_NAME
    );

    /**
     * The conf threshold above which variants are not included in the position tables.
     * This value is used to construct the genotype information of those missing samples
     * when they are merged together into a {@link VariantContext} object
     */
    public static final int MISSING_CONF_THRESHOLD = 60;

    private final boolean precomputedResultsMode;

    private final VariantContextWriter vcfWriter;

    private final VCFHeader vcfHeader;

    private final ReferenceDataSource refSource;

    private final ReferenceConfidenceVariantContextMerger variantContextMerger;

    private final GnarlyGenotyperEngine gnarlyGenotyper;

    private final VariantAnnotatorEngine variantAnnotator;

    private final String projectID;

    /** Set of sample names seen in the variant data from BigQuery. */
    private final Set<String> sampleNames = new HashSet<>();

    /**
     * Map between contig name and the BigQuery table containing position data from that contig.
     */
    private final Map<String, String> contigToPositionTableMap;

    /**
     * Map between contig name and the BigQuery table containing variant data from that contig.
     */
    private final Map<String, String> contigToVariantTableMap;

    /**
     * Map between contig name and the BigQuery table containing allele data from that contig.
     */
    private final Map<String, String> contigToAltAlleleTableMap;

    /**
     * Map between contig name and the BigQuery table containing the list of samples
     */
    private final Map<String, String> contigToSampleTableMap;
    
    private final int queryRecordLimit;

    private final boolean useOptimizedQuery;
    private final boolean useCohortExtractQuery;
    private final String filteringFQTableName;
    private final boolean useModelFeatureExtractQuery;
    private final boolean trainingSitesOnly;

    private final boolean doLocalSort;

    private final int localSortMaxRecordsInRam;

    private final boolean runQueryOnly;

    private final boolean disableGnarlyGenotyper;
    private final boolean enableVariantAnnotator;

    private final boolean runQueryInBatchMode;

    private final boolean printDebugInformation;

    private double vqsLodSNPThreshold = 0.0;
    private double vqsLodINDELThreshold = 0.0;

    private final ProgressMeter progressMeter;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

    EvoquerEngine( final String projectID,
                   final Map<String, Evoquer.EvoquerDataset> datasetMap,
                   final int queryRecordLimit,
                   final boolean useOptimizedQuery,
                   final VariantContextWriter vcfWriter,
                   final Set<VCFHeaderLine> toolDefaultVCFHeaderLines,
                   final VariantAnnotatorEngine annotationEngine,
                   final ReferenceDataSource refSource,
                   final String sampleTableName,
                   final boolean doLocalSort,
                   final boolean useCohortExtractQuery,
                   final String filteringFQTableName,
                   final boolean useModelFeatureExtractQuery,
                   final boolean  trainingSitesOnly,
                   final int localSortMaxRecordsInRam,
                   final boolean runQueryOnly,
                   final boolean disableGnarlyGenotyper,
                   final boolean enableVariantAnnotator,
                   final boolean keepAllSitesInGnarlyGenotyper,
                   final boolean runQueryInBatchMode,
                   final boolean printDebugInformation,
                   final double vqsLodSNPThreshold,
                   final double vqsLodINDELThreshold,
                   final ProgressMeter progressMeter ) {

        // We were given a dataset map, so we're going to do live queries against BigQuery
        this.precomputedResultsMode = false;
        this.useOptimizedQuery = useOptimizedQuery;
        this.useCohortExtractQuery = useCohortExtractQuery;
        this.filteringFQTableName = filteringFQTableName;
        this.useModelFeatureExtractQuery = useModelFeatureExtractQuery;
        this.trainingSitesOnly = trainingSitesOnly;

        this.doLocalSort = doLocalSort;
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.projectID = projectID;
        this.queryRecordLimit = queryRecordLimit;
        this.runQueryOnly = runQueryOnly;
        this.disableGnarlyGenotyper = disableGnarlyGenotyper;
        this.enableVariantAnnotator = enableVariantAnnotator;
        this.runQueryInBatchMode = runQueryInBatchMode;
        this.printDebugInformation = printDebugInformation;
        this.vqsLodSNPThreshold = vqsLodSNPThreshold;
        this.vqsLodINDELThreshold = vqsLodINDELThreshold;
        this.progressMeter = progressMeter;

        final Map<String, String> tmpContigToPositionTableMap = new HashMap<>();
        final Map<String, String> tmpContigToVariantTableMap = new HashMap<>();
        final Map<String, String> tmpContigToAltAlleleTableMap = new HashMap<>();
        final Map<String, String> tmpContigToSampleTableMap = new HashMap<>();
        for ( final Map.Entry<String, Evoquer.EvoquerDataset> datasetEntry : datasetMap.entrySet() ) {
            tmpContigToPositionTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + datasetEntry.getValue().getPetTableName());
            tmpContigToVariantTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + datasetEntry.getValue().getVetTableName());
            tmpContigToAltAlleleTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + datasetEntry.getValue().getAltAlleleTableName());
            tmpContigToSampleTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + sampleTableName);
        }
        contigToPositionTableMap = Collections.unmodifiableMap(tmpContigToPositionTableMap);
        contigToVariantTableMap = Collections.unmodifiableMap(tmpContigToVariantTableMap);
        contigToAltAlleleTableMap = Collections.unmodifiableMap(tmpContigToAltAlleleTableMap);
        contigToSampleTableMap = Collections.unmodifiableMap(tmpContigToSampleTableMap);

        // Get the samples used in the dataset:
        if (!useModelFeatureExtractQuery) {
            populateSampleNames();
        }
        this.vcfHeader = generateVcfHeader(toolDefaultVCFHeaderLines, refSource.getSequenceDictionary());

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);
        this.gnarlyGenotyper = new GnarlyGenotyperEngine(keepAllSitesInGnarlyGenotyper, GenotypeCalculationArgumentCollection.DEFAULT_MAX_ALTERNATE_ALLELES, false, false);
        this.variantAnnotator = new VariantAnnotatorEngine(Collections.singletonList(new ChromosomeCounts()), null, Collections.emptyList(), false, false);

    }

    EvoquerEngine( final List<String> sampleNames,
                   final VariantContextWriter vcfWriter,
                   final Set<VCFHeaderLine> toolDefaultVCFHeaderLines,
                   final VariantAnnotatorEngine annotationEngine,
                   final ReferenceDataSource refSource,
                   final boolean disableGnarlyGenotyper,
                   final boolean keepAllSitesInGnarlyGenotyper,
                   final boolean printDebugInformation,
                   final double vqsLodSNPThreshold,
                   final double vqsLodINDELThreshold,
                   final ProgressMeter progressMeter ) {

        // We weren't given a dataset map, so we're not going to do any live queries against BigQuery,
        // and will instead expect to be given URIs to precomputed results.
        this.precomputedResultsMode = true;
        this.useOptimizedQuery = false;
        this.useCohortExtractQuery = false;
        this.filteringFQTableName = null;
        this.useModelFeatureExtractQuery = false;
        this.trainingSitesOnly = false;

        this.doLocalSort = false;
        this.localSortMaxRecordsInRam = 0;

        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.disableGnarlyGenotyper = disableGnarlyGenotyper;
        this.enableVariantAnnotator = false;
        this.printDebugInformation = printDebugInformation;
        this.vqsLodSNPThreshold = vqsLodSNPThreshold;
        this.vqsLodINDELThreshold = vqsLodINDELThreshold;
        this.progressMeter = progressMeter;

        this.projectID = null;
        this.contigToPositionTableMap = null;
        this.contigToVariantTableMap = null;
        this.contigToAltAlleleTableMap = null;
        this.contigToSampleTableMap = null;
        this.queryRecordLimit = 0;
        this.runQueryOnly = false;
        this.runQueryInBatchMode = false;

        this.sampleNames.addAll(sampleNames);
        this.vcfHeader = generateVcfHeader(toolDefaultVCFHeaderLines, refSource.getSequenceDictionary());
        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);
        this.gnarlyGenotyper = new GnarlyGenotyperEngine(keepAllSitesInGnarlyGenotyper, GenotypeCalculationArgumentCollection.DEFAULT_MAX_ALTERNATE_ALLELES, false, false);
        this.variantAnnotator = new VariantAnnotatorEngine(Collections.singletonList(new ChromosomeCounts()), null, Collections.emptyList(), false, false);

    }

    /**
     * Connects to the BigQuery table for the given interval and pulls out the information on the samples that
     * contain variants.
     * @param interval {@link SimpleInterval} over which to query the BigQuery table.
     */
    void evokeInterval(final SimpleInterval interval) {
        if ( precomputedResultsMode ) {
            throw new GATKException("Cannot do live queries in precomputedResultsMode");
        }

        if (!contigToPositionTableMap.containsKey(interval.getContig())) {
            logger.warn("Contig missing from contigPositionExpandedTableMap, ignoring interval: " + interval.toString());
            return;
        }

        if ( useCohortExtractQuery ) {
            if ( filteringFQTableName == null || filteringFQTableName.equals("") ) {
                logger.warn("--cohort-extract-filter-table must be specified when extracting a cohort! ");
                return;
            }

            final String variantQueryString = getCohortExtractQueryString(interval);
            final List<String> fieldsToRetrieve = Arrays.asList("position", "sample", "state", "ref", "alt", "call_YNG_STATUS", "call_AS_VQSLOD", "call_GT", "call_GQ", "call_RGQ");
            final BigQueryUtils.StorageAPIAvroReader storageAPIAvroReader = BigQueryUtils.executeQueryWithStorageAPI(variantQueryString, fieldsToRetrieve, projectID, runQueryInBatchMode);

            createVariantsFromUngroupedTableResult(storageAPIAvroReader, interval.getContig());
        }  else if ( useModelFeatureExtractQuery ) {
            final String featureQueryString  = getVQSRFeatureExtractQueryString(interval, trainingSitesOnly);
            logger.info(featureQueryString);
            final List<String> fieldsToRetrieve = Arrays.asList("position", "ref", "allele", "RAW_QUAL", "ref_ad", "AS_MQRankSum", "AS_MQRankSum_ft", "AS_ReadPosRankSum", "AS_ReadPosRankSum_ft", "RAW_MQ", "RAW_AD", "RAW_AD_GT_1", "SB_REF_PLUS","SB_REF_MINUS","SB_ALT_PLUS","SB_ALT_MINUS");
            final BigQueryUtils.StorageAPIAvroReader storageAPIAvroReader = BigQueryUtils.executeQueryWithStorageAPI(featureQueryString, fieldsToRetrieve, projectID, runQueryInBatchMode);

            createVQSRInputFromTableResult(storageAPIAvroReader, interval.getContig());
        } else {
            String variantQueryString;
            List<String> fieldsToRetrieve;

            if ( doLocalSort ) {
                variantQueryString = useOptimizedQuery ?
                        getOptimizedUngroupedVariantQueryString(interval) :
                        getUngroupedVariantQueryString(interval);
                // TODO: get the selected fields from the query string itself somehow
                fieldsToRetrieve = Arrays.asList("position", "sample", "state", "ref", "alt", "AS_RAW_MQ", "AS_RAW_MQRankSum", "AS_QUALapprox", "AS_RAW_ReadPosRankSum", "AS_SB_TABLE", "AS_VarDP", "call_GT", "call_AD", "call_DP", "call_GQ", "call_PGT", "call_PID", "call_PL");
            } else {
                variantQueryString = useOptimizedQuery ?
                        getOptimizedGroupedVariantQueryString(interval) :
                        getGroupedVariantQueryString(interval);
                // TODO: get the selected fields from the query string itself somehow
                fieldsToRetrieve = Arrays.asList("position", "values");
            }

            // Execute the query:
            final BigQueryUtils.StorageAPIAvroReader storageAPIAvroReader = BigQueryUtils.executeQueryWithStorageAPI(variantQueryString, fieldsToRetrieve, projectID, runQueryInBatchMode);

            if ( doLocalSort ) {
                createVariantsFromUngroupedTableResult(storageAPIAvroReader, interval.getContig());
            } else {
                createVariantsFromGroupedTableResult(storageAPIAvroReader, interval.getContig());
            }
        }
    }

    void evokeAvroResult(final String avroResultFile, final String contig) {
        if ( ! precomputedResultsMode ) {
            throw new GATKException("Must be in precomputed results mode to accept precomputed avro inputs");
        }

        final GATKAvroReader avroReader = new GCSAvroReader(avroResultFile);
        createVariantsFromGroupedTableResult(avroReader, contig);
    }

    /**
     * Generates a {@link VCFHeader} object based on the VariantContext objects queried from the BigQuery backend.
     * If no objects have been queried, this will return a default {@link VCFHeader}.
     * @param defaultHeaderLines The default header lines to be added to the top of the VCF header.
     * @param sequenceDictionary The SequenceDictionary of the reference on which the variants are based.
     * @return A {@link VCFHeader} object representing the header for all variants that have been queried from the BigQuery backend.
     */
    private VCFHeader generateVcfHeader(final Set<VCFHeaderLine> defaultHeaderLines,
                                       final SAMSequenceDictionary sequenceDictionary) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        headerLines.addAll( getEvoquerVcfHeaderLines() );
        headerLines.addAll( defaultHeaderLines );

        final VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(sequenceDictionary);

        return header;
    }

    VCFHeader getHeader() {
        return vcfHeader;
    }

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }

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
                VCFConstants.ALLELE_NUMBER_KEY,
                VCFConstants.END_KEY
        );

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS
        );

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
        
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_SB_TABLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VQS_LOD_KEY));

        // TODO: Temporary.  We don't really want these as FORMAT fields,
        headerLines.add(
                new VCFInfoHeaderLine(AS_VQS_LOD_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each alt allele, the log odds of being a true variant versus being false under the trained gaussian mixture model")
        );
        headerLines.add(
                new VCFInfoHeaderLine(AS_YNG_STATUS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each alt allele, status of the YNG filter")
        );

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_FISHER_STRAND_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.FISHER_STRAND_KEY));


        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXCESS_HET_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));

        // TODO: There m ust be a more appropriate constant to use for these
        headerLines.add(new VCFFilterHeaderLine(PASSES_FILTERS_v4, "PASSING"));
        headerLines.add(new VCFFilterHeaderLine("NAY", "Site is Nay in the YNG table"));
        headerLines.add(new VCFFilterHeaderLine("VQSRTranchSNP", "Site fails to exceed the SNP tranch threshold"));
        headerLines.add(new VCFFilterHeaderLine("VQSRTrachINDEL", "Site fails to exceel the INDEL tranch threshold"));

        return headerLines;
    }

    private String getPositionTableForContig(final String contig ) {
        return contigToPositionTableMap.get(contig);
    }

    private String getVariantTableForContig(final String contig ) {
        return contigToVariantTableMap.get(contig);
    }

    private String getTableQualifier() {
        return projectID;
    }

    private String getFQTableName( final String tableName ) {
        return getTableQualifier() + "." + tableName;
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the position
     * data specified in the given {@code interval}.
     *
     * Uses {@link #projectID} for the project of the BigQuery table.
     * Assumes the tables have dataset information in them.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private String getFQPositionTable(final SimpleInterval interval) {
        return getFQTableName(getPositionTableForContig( interval.getContig() ));
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the variant
     * data specified in the given {@code interval}.
     *
     * Uses {@link #projectID} for the project of the BigQuery table.
     * Assumes the tables have dataset information in them.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private String getFQVariantTable(final SimpleInterval interval) {
        return getFQTableName(getVariantTableForContig( interval.getContig() ));
    }

    private void createVariantsFromGroupedTableResult( final GATKAvroReader avroReader, final String contig ) {
        final org.apache.avro.Schema schema = avroReader.getSchema();

        final Set<String> columnNames = new HashSet<>();
        if ( schema.getField(POSITION_FIELD_NAME) == null || schema.getField(VALUES_ARRAY_FIELD_NAME) == null ) {
            throw new UserException("Records must contain position and values columns");
        }
        schema.getField(VALUES_ARRAY_FIELD_NAME).schema().getElementType().getFields().forEach(field -> columnNames.add(field.name()));
        validateSchema(columnNames);

        for ( final GenericRecord row : avroReader ) {
            if ( runQueryOnly ) {
                continue;
            }

            ++totalNumberOfSites;

            final long currentPosition = Long.parseLong(row.get(POSITION_FIELD_NAME).toString());

            if ( printDebugInformation ) {
                logger.info(contig + ":" + currentPosition + ": found record: " + row);
            }

            @SuppressWarnings("unchecked")
            final GenericData.Array<GenericRecord> sampleRecords = (GenericData.Array<GenericRecord>) row.get(VALUES_ARRAY_FIELD_NAME);

            processSampleRecordsForPosition(currentPosition, contig, sampleRecords, columnNames);
        }
    }

    private EvoquerSortingCollection<GenericRecord> getEvoquerSortingCollection(org.apache.avro.Schema schema) {
        final EvoquerSortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        final Comparator<GenericRecord> sortingCollectionComparator = new Comparator<GenericRecord>() {
            @Override
            public int compare( GenericRecord o1, GenericRecord o2 ) {
                final long firstPosition = Long.parseLong(o1.get(POSITION_FIELD_NAME).toString());
                final long secondPosition = Long.parseLong(o2.get(POSITION_FIELD_NAME).toString());

                return Long.compare(firstPosition, secondPosition);
            }
        };
        return EvoquerSortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, sortingCollectionComparator, localSortMaxRecordsInRam);
    }

    private void createVQSRInputFromTableResult(final GATKAvroReader avroReader, final String contig) {
        final org.apache.avro.Schema schema = avroReader.getSchema();
        final Set<String> columnNames = new HashSet<>();
        if ( schema.getField(POSITION_FIELD_NAME) == null ) {
            throw new UserException("Records must contain a position column");
        }
        schema.getFields().forEach(field -> columnNames.add(field.name()));

        // TODO: this hardcodes a list of required fields... which this case doesn't require!
        // should genericize/parameterize so we can also validate the schema here
        // validateSchema(columnNames);

        EvoquerSortingCollection<GenericRecord> sortingCollection =  getEvoquerSortingCollection(schema);

        for ( final GenericRecord queryRow : avroReader ) {
            if ( runQueryOnly ) {
                continue;
            }

            sortingCollection.add(queryRow);
        }

        sortingCollection.printTempFileStats();

        if ( runQueryOnly ) {
            return;
        }

        for ( final GenericRecord row : sortingCollection ) {
            processVQSRRecordForPosition(row, contig);
        }
    }

    private void processVQSRRecordForPosition(GenericRecord rec, String contig) {
        final long position = Long.parseLong(rec.get(POSITION_FIELD_NAME).toString());

        // I don't understand why the other modes iterate through a list of all columns, and then
        //  switch based on the column name rather than just getting the columns desired directly (like I'm going to do here).
        // It might be because those field names are just prefixed versions of standard VCF fields?

        // TODO: de-python names (no  _)
        String ref = rec.get("ref").toString();
        String allele = rec.get("allele").toString();

        // Numbers are returned as Long (sci notation)
        Double qual = Double.valueOf(rec.get("RAW_QUAL").toString());

        Object o_raw_ref_ad = rec.get("ref_ad");
        Double raw_ref_ad = (o_raw_ref_ad==null)?0:Double.valueOf(o_raw_ref_ad.toString());

        Object o_AS_MQRankSum = rec.get("AS_MQRankSum");
        Float AS_MQRankSum = (o_AS_MQRankSum==null)?null:Float.parseFloat(o_AS_MQRankSum.toString());

        Object o_AS_ReadPosRankSum = rec.get("AS_ReadPosRankSum");
        Float AS_ReadPosRankSum = (o_AS_ReadPosRankSum==null)?null:Float.parseFloat(o_AS_ReadPosRankSum.toString());

        Double raw_mq = Double.valueOf(rec.get("RAW_MQ").toString());
        Double raw_ad = Double.valueOf(rec.get("RAW_AD").toString());
        Double raw_ad_gt_1 = Double.valueOf(rec.get("RAW_AD_GT_1").toString());

        // TODO: KCIBUL QUESTION -- if we skip this... we won't have YNG Info @ extraction time?
        if (raw_ad == 0) {
            logger.info("skipping " + position + " because it has no alternate reads!");
            return;
        }

        int sb_ref_plus = Double.valueOf(rec.get("SB_REF_PLUS").toString()).intValue();
        int sb_ref_minus = Double.valueOf(rec.get("SB_REF_MINUS").toString()).intValue();
        int sb_alt_plus = Double.valueOf(rec.get("SB_ALT_PLUS").toString()).intValue();
        int sb_alt_minus = Double.valueOf(rec.get("SB_ALT_MINUS").toString()).intValue();


//        logger.info("processing " + contig + ":" + position);



        // NOTE: if VQSR required a merged VCF (e.g. multiple alleles on a  given row) we have to do some merging here...

        final VariantContextBuilder builder = new VariantContextBuilder();

        builder.chr(contig);
        builder.start(position);

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(allele, false));
        builder.alleles(alleles);
        builder.stop(position + alleles.get(0).length() - 1);


        double qd_depth = (raw_ref_ad + raw_ad_gt_1);
        double as_qd = QualByDepth.fixTooHighQD( (qual /  qd_depth) );

        final int[][] refAltTable = new int[][]{new int[]{sb_ref_plus, sb_ref_minus}, new int[]{ sb_alt_plus, sb_alt_minus}};

        double fs = QualityUtils.phredScaleErrorRate(Math.max(FisherStrand.pValueForContingencyTable(refAltTable), AS_StrandBiasTest.MIN_PVALUE));
        double sor = StrandOddsRatio.calculateSOR(refAltTable);

        double mq = Math.sqrt( raw_mq / raw_ad);

        builder.attribute("AS_QD", String.format("%.2f", as_qd) );
        builder.attribute("AS_FS", String.format("%.3f", fs));
        builder.attribute("AS_MQ", String.format("%.2f", mq) );
        builder.attribute("AS_MQRankSum", AS_MQRankSum==null?".":String.format("%.3f", AS_MQRankSum) );
        builder.attribute("AS_ReadPosRankSum", AS_ReadPosRankSum==null?".":String.format("%.3f", AS_ReadPosRankSum));
        builder.attribute("AS_SOR", String.format("%.3f", sor));

        // check out 478765 -- we need to "merge" different variant  contexts  at the same position.
        // I think if we just set the right "base" annotations it's possible the standard annotation processing
        // merging will take care of it?
        // Also -- is it possible to write a JS UDF like areEqual(ref1, alt2, ref2, alt2)?
        VariantContext vc = builder.make();
        vcfWriter.add(vc);
        progressMeter.update(vc);


    }

    private void createVariantsFromUngroupedTableResult(final GATKAvroReader avroReader, final String contig) {

        final org.apache.avro.Schema schema = avroReader.getSchema();

        final Set<String> columnNames = new HashSet<>();
        if ( schema.getField(POSITION_FIELD_NAME) == null ) {
            throw new UserException("Records must contain a position column");
        }
        schema.getFields().forEach(field -> columnNames.add(field.name()));
        validateSchema(columnNames);

        EvoquerSortingCollection<GenericRecord> sortingCollection =  getEvoquerSortingCollection(schema);

        for ( final GenericRecord queryRow : avroReader ) {
            if ( runQueryOnly ) {
                continue;
            }

            sortingCollection.add(queryRow);
        }

        sortingCollection.printTempFileStats();

        if ( runQueryOnly ) {
            return;
        }

        final List<GenericRecord> currentPositionRecords = new ArrayList<>(sampleNames.size() * 2);
        long currentPosition = -1;

        for ( final GenericRecord sortedRow : sortingCollection ) {
            final long rowPosition = Long.parseLong(sortedRow.get(POSITION_FIELD_NAME).toString());

            if ( rowPosition != currentPosition && currentPosition != -1 ) {
                ++totalNumberOfSites;
                processSampleRecordsForPosition(currentPosition, contig, currentPositionRecords, columnNames);

                currentPositionRecords.clear();
            }

            currentPositionRecords.add(sortedRow);
            currentPosition = rowPosition;
        }

        if ( ! currentPositionRecords.isEmpty() ) {
            ++totalNumberOfSites;
            processSampleRecordsForPosition(currentPosition, contig, currentPositionRecords, columnNames);
        }
    }

    private void processSampleRecordsForPosition(final long currentPosition, final String contig, final Iterable<GenericRecord> sampleRecordsAtPosition, final Set<String> columnNames) {
        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final Set<String> currentPositionSamplesSeen = new HashSet<>();
        boolean currentPositionHasVariant = false;
        final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, currentPosition, currentPosition).getBaseString(), true);
        int numRecordsAtPosition = 0;

        final HashMap<Allele, HashMap<Allele, Double>> vqsLodMap = new HashMap<>();
        final HashMap<Allele, HashMap<Allele, String>> yngMap = new HashMap<>();

        for ( final GenericRecord sampleRecord : sampleRecordsAtPosition ) {
            final String sampleName = sampleRecord.get(SAMPLE_FIELD_NAME).toString();
            currentPositionSamplesSeen.add(sampleName);
            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + currentPosition + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            switch (sampleRecord.get(STATE_FIELD_NAME).toString()) {
                case "v":   // Variant
                    ++totalNumberOfVariants;
                    unmergedCalls.add(createVariantContextFromSampleRecord(sampleRecord, columnNames, contig, currentPosition, sampleName, vqsLodMap, yngMap));
                    currentPositionHasVariant = true;
                    break;
                case "0":   // Non Variant Block with GQ < 10
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 0));
                    break;
                case "1":  // Non Variant Block with 10 <=  GQ < 20
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 10));
                    break;
                case "2":  // Non Variant Block with 20 <= GQ < 30
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 20));
                    break;
                case "3":  // Non Variant Block with 30 <= GQ < 40
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 30));
                    break;
                case "4":  // Non Variant Block with 40 <= GQ < 50
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 40));
                    break;
                case "5":  // Non Variant Block with 50 <= GQ < 60
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 50));
                    break;
                case "6":  // Non Variant Block with 60 <= GQ (usually omitted from tables)
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 60));
                    break;
                case "*":   // Spanning Deletion: TODO handle correctly
                    break;
                case "m":   // Missing
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                default:
                    throw new GATKException("Unrecognized state: " + sampleRecord.get(STATE_FIELD_NAME).toString());
            }

        }

        if ( printDebugInformation ) {
            logger.info(contig + ":" + currentPosition + ": processed " + numRecordsAtPosition + " total sample records");
        }

        finalizeCurrentVariant(unmergedCalls, currentPositionSamplesSeen, currentPositionHasVariant, contig, currentPosition, refAllele, vqsLodMap, yngMap);
    }

    private void validateSchema(final Set<String> columnNames) {
        for ( final String requiredField : REQUIRED_FIELDS ) {
            if ( ! columnNames.contains(requiredField) ) {
                throw new UserException("Missing required column: " + requiredField +
                    ". Actual columns encountered were: " + columnNames);
            }
        }
    }

    // vqsLogMap and yngMap are in/out parameters for this method. i.e. they are modified by this method
    private VariantContext createVariantContextFromSampleRecord(final GenericRecord sampleRecord, final Set<String> columnNames, final String contig, final long startPosition, final String sample, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(startPosition);

        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(sampleRecord.get(REF_ALLELE_FIELD_NAME).toString(), true);
        alleles.add(ref);
        List<Allele> altAlleles = Arrays.stream(sampleRecord.get(ALT_ALLELE_FIELD_NAME).toString().split(MULTIVALUE_FIELD_DELIMITER))
                .map(altAllele -> Allele.create(altAllele, false)).collect(Collectors.toList());
        alleles.addAll(altAlleles);
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sample);

        vqsLodMap.putIfAbsent(ref, new HashMap<>());
        yngMap.putIfAbsent(ref, new HashMap<>());


        for ( final String columnName : columnNames ) {
            if ( REQUIRED_FIELDS.contains(columnName) || columnName.equals(POSITION_FIELD_NAME) ) {
                continue;
            }

            final Object columnValue = sampleRecord.get(columnName);
            if ( columnValue == null ) {
                continue;
            }
            final String columnValueString = columnValue.toString();

            if ( columnName.startsWith(GENOTYPE_FIELD_PREFIX) ) {
                final String genotypeAttributeName = columnName.substring(GENOTYPE_FIELD_PREFIX.length());

                if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_KEY) ) {
                    final List<Allele> genotypeAlleles =
                            Arrays.stream(columnValueString.split("[/|]"))
                            .map(Integer::parseInt)
                            .map(alleleIndex -> alleles.get(alleleIndex))
                            .collect(Collectors.toList());
                    genotypeBuilder.alleles(genotypeAlleles);
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    genotypeBuilder.GQ(Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_PL_KEY) ) {
                    genotypeBuilder.PL(Arrays.stream(columnValueString.split(MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(VCFConstants.DEPTH_KEY) ) {
                    genotypeBuilder.DP(Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY) ) {
                    genotypeBuilder.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS) ) {
                    genotypeBuilder.AD(Arrays.stream(columnValueString.split(MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(GATKVCFConstants.AS_VQS_LOD_KEY) ) {
                    HashMap<Allele, Double> innerVqslodMap = vqsLodMap.get(ref);
                    double[] vqslodValues = Arrays.stream(columnValueString.split(MULTIVALUE_FIELD_DELIMITER)).mapToDouble(Double::parseDouble).toArray();
                    // should we use put or putIfAbsent - this may insert the same value multiple times. same with YNG below
                    new IndexRange(0, vqslodValues.length).forEach(i -> innerVqslodMap.put(altAlleles.get(i), vqslodValues[i]));
                } else if ( genotypeAttributeName.equals("YNG_STATUS") ) {
                    HashMap<Allele, String> innerYNGMap = yngMap.get(ref);
                    String[] yngValues = columnValueString.split(MULTIVALUE_FIELD_DELIMITER);
                    new IndexRange(0, yngValues.length).forEach(i -> innerYNGMap.put(altAlleles.get(i), yngValues[i]));
                } else {
                    genotypeBuilder.attribute(genotypeAttributeName, columnValueString);
                }
            } else {
                builder.attribute(columnName, columnValueString);
            }
        }

        builder.genotypes(genotypeBuilder.make());

        return builder.make();
    }

    private VariantContext createRefSiteVariantContext(final String sample, final String contig, final long start, final Allele refAllele, final int gq) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(start);

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(refAllele);
        alleles.add(Allele.NON_REF_ALLELE);
        builder.alleles(alleles);

        builder.stop(start);

        genotypeBuilder.name(sample);

        builder.attribute(VCFConstants.END_KEY, Long.toString(start));

        final List<Allele> genotypeAlleles = new ArrayList<>();
        genotypeAlleles.add(refAllele);
        genotypeAlleles.add(refAllele);
        genotypeBuilder.alleles(genotypeAlleles);

        genotypeBuilder.GQ(gq);

        builder.genotypes(genotypeBuilder.make());
        
        return builder.make();
    }

    private void finalizeCurrentVariant(final List<VariantContext> unmergedCalls, final Set<String> currentVariantSamplesSeen, final boolean currentPositionHasVariant, final String contig, final long start, final Allele refAllele, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap) {
        // If there were no variants at this site, we don't emit a record and there's nothing to do here
        if ( ! currentPositionHasVariant ) {
            return;
        }

        // Find missing samples and synthesize GQ 60s
        final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        for ( final String missingSample : samplesNotEncountered ) {
            unmergedCalls.add(createRefSiteVariantContext(missingSample, contig, start, refAllele, MISSING_CONF_THRESHOLD));
        }

        // Note that we remove NON_REF in our variantContextMerger if the GnarlyGenotyper is disabled, but
        // keep NON_REF around if the GnarlyGenotyper is enabled. This is because the GnarlyGenotyper expects
        // the NON_REF allele to still be present, and will give incorrect results if it's not.
        final VariantContext mergedVC = variantContextMerger.merge(unmergedCalls, new SimpleInterval(contig, (int) start, (int) start), refAllele.getBases()[0], disableGnarlyGenotyper, false, true);


        LinkedHashMap<Allele, Double> remappedVqsLodMap = remapAllelesInMap(mergedVC, vqsLodMap);
        LinkedHashMap<Allele, String> remappedYngMap = remapAllelesInMap(mergedVC, yngMap);

        final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);
        builder.attribute(GATKVCFConstants.AS_VQS_LOD_KEY, remappedVqsLodMap.values() );
        builder.attribute(GATKVCFConstants.AS_YNG_STATUS_KEY, remappedYngMap.values() );

        int refLength = mergedVC.getReference().length();

        // TODO: should we use the FilterVariantTranches tool instead of reimplementing it ourselves???

        // if there are any Yays, the site is PASS
        if (remappedYngMap.values().contains("Y")) {
            builder.filter("PASS");
        } else if (remappedYngMap.values().contains("N")) {
            // TODO: do we want to remove this variant?
              builder.filter("NAY");
        } else {
            if (remappedYngMap.values().contains("G")) {
                Optional<Double> snpMax = remappedVqsLodMap.entrySet().stream().filter(entry -> entry.getKey().length() == refLength).map(entry -> entry.getValue()).max(Double::compareTo);
                if (snpMax.isPresent() && snpMax.get() < vqsLodSNPThreshold) {
                    // TODO: add in sensitivities
                    builder.filter("VQSRTranchSNP");
                }
                Optional<Double> indelMax = remappedVqsLodMap.entrySet().stream().filter(entry -> entry.getKey().length() != refLength).map(entry -> entry.getValue()).max(Double::compareTo);
                if (indelMax.isPresent() && indelMax.get() < vqsLodINDELThreshold) {
                    // TODO: add in sensitivities
                    builder.filter("VQSRTrancheINDEL");
                }
            }
            // TODO: what if there is nothing in the YNG table?
        }

        final VariantContext filteredVC = builder.make();

        final VariantContext finalizedVC = disableGnarlyGenotyper ? filteredVC : gnarlyGenotyper.finalizeGenotype(filteredVC);
        final VariantContext annotatedVC = enableVariantAnnotator ?
                variantAnnotator.annotateContext(finalizedVC, new FeatureContext(), null, null, a -> true): finalizedVC;

        if ( annotatedVC != null ) {
            vcfWriter.add(annotatedVC);
            progressMeter.update(annotatedVC);
        }
        else if ( finalizedVC != null ) { // GnarlyGenotyper returns null for variants it refuses to output
            vcfWriter.add(finalizedVC);
            progressMeter.update(finalizedVC);
        }
        else {
            logger.warn(String.format("GnarlyGenotyper returned null for site %s:%s", contig, start));
            progressMeter.update(mergedVC);
        }
    }

    private <T> LinkedHashMap<Allele, T> remapAllelesInMap(VariantContext vc, HashMap<Allele, HashMap<Allele, T>> datamap) {
        // get the extended reference
        Allele ref = vc.getReference();

        // create ordered results map
        LinkedHashMap<Allele, T> results = new LinkedHashMap<>();
        vc.getAlternateAlleles().stream().forEachOrdered(allele -> results.put(allele, null));

        List<Allele> newAlleles = new ArrayList<>();
        datamap.entrySet().stream().forEachOrdered(entry -> {
            if (entry.getKey() == ref) {
                // reorder
                entry.getValue().entrySet().stream().forEach(altMapEntry -> results.put(altMapEntry.getKey(), altMapEntry.getValue()));
            } else {
                // remap
                List<Allele> allAlleles = new ArrayList<>();
                allAlleles.add(entry.getKey());
                allAlleles.addAll(entry.getValue().keySet());
                VariantContextBuilder vcb = new VariantContextBuilder(vc.getSource(), vc.getContig(), vc.getStart(), vc.getEnd(), allAlleles);
                Map<Allele, Allele> alleleMapping = GATKVariantContextUtils.createAlleleMapping(ref, vcb.make(), newAlleles);
                alleleMapping.entrySet().stream().forEach(mapped -> results.put(mapped.getValue(), entry.getValue().get(mapped.getKey())));
            }
        });
        return results;
    }

    // HACK to deal with malformed allele-specific annotations in our test dataset
    private boolean alleleSpecificAnnotationsAreMalformed( final VariantContext vc ) {
        // AS_RAW_MQ  AS_RAW_MQRankSum  AS_QUALapprox  AS_RAW_ReadPosRankSum  AS_SB_TABLE  AS_VarDP
        final int AS_RAW_MQ_values = StringUtils.countMatches(vc.getAttributeAsString("AS_RAW_MQ", ""), "|") + 1;
        final int AS_RAW_MQRankSum_values = StringUtils.countMatches(vc.getAttributeAsString("AS_RAW_MQRankSum", ""), "|") + 1;
        final int AS_QUALapprox_values = StringUtils.countMatches(vc.getAttributeAsString("AS_QUALapprox", ""), "|") + 1;
        final int AS_RAW_ReadPosRankSum_values = StringUtils.countMatches(vc.getAttributeAsString("AS_RAW_ReadPosRankSum", ""), "|") + 1;
        final int AS_SB_TABLE_values = StringUtils.countMatches(vc.getAttributeAsString("AS_SB_TABLE", ""), "|") + 1;
        final int AS_VarDP_values = StringUtils.countMatches(vc.getAttributeAsString("AS_VarDP", ""), "|") + 1;

        if ( AS_RAW_MQ_values != vc.getNAlleles() ||
             AS_RAW_MQRankSum_values != vc.getNAlleles() ||
             AS_QUALapprox_values != vc.getNAlleles() - 1 ||
             AS_RAW_ReadPosRankSum_values != vc.getNAlleles() ||
             AS_SB_TABLE_values != vc.getNAlleles() ||
             AS_VarDP_values != vc.getNAlleles() ) {

            return true;
        }

        return false;
    }

    private void populateSampleNames() {
        // TODO: For now, use the sample list from an arbitrary dataset's sample_list table.
        // TODO: Eventually, may need to crosscheck sample lists across datasets
        final String sampleTableName = contigToSampleTableMap.entrySet().iterator().next().getValue();

        // Get the query string:
        final String sampleListQueryString = getSampleListQueryString(sampleTableName);
        
        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        // Show our pretty results:
        if ( printDebugInformation ) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info("\n" + prettyQueryResults);
        }

        // Add our samples to our map:
        for ( final FieldValueList row : result.iterateAll() ) {
            sampleNames.add( row.get(0).getStringValue() );
        }
    }

    private String getGroupedVariantQueryString( final SimpleInterval interval ) {
        String limitString = "";
        if ( queryRecordLimit > 0 ) {
            limitString = "LIMIT " + queryRecordLimit;
        }
        
        return String.format(
                "SELECT position, ARRAY_AGG(STRUCT( sample, state, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, call_GT, call_AD, call_DP, call_GQ, call_PGT, call_PID, call_PL  )) AS values\n" +
                "FROM `%s` AS pet\n" +
                "LEFT OUTER JOIN `%s` AS vet\n" +
                "USING (position, sample)\n" +
                "WHERE (position >= %d AND position <= %d)\n" +
                "GROUP BY position\n" +
                limitString,
                getFQPositionTable(interval),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd());
    }

    private String getUngroupedVariantQueryString( final SimpleInterval interval ) {
        String limitString = "";
        if ( queryRecordLimit > 0 ) {
            limitString = "LIMIT " + queryRecordLimit;
        }

        return String.format(
                "SELECT position, sample, state, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, call_GT, call_AD, call_DP, call_GQ, call_PGT, call_PID, call_PL\n" +
                        "FROM `%s` AS pet\n" +
                        "LEFT OUTER JOIN `%s` AS vet\n" +
                        "USING (position, sample)\n" +
                        "WHERE (position >= %d AND position <= %d)\n" +
                        limitString,
                getFQPositionTable(interval),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd());
    }

    private String getOptimizedGroupedVariantQueryString(final SimpleInterval interval) {
        String limitString = "";
        if (queryRecordLimit > 0) {
            limitString = "LIMIT " + queryRecordLimit;
        }

        return String.format(
                "WITH new_pet AS\n" +
                        "(\n" +
                        "  SELECT * FROM `%s`\n" +
                        "  WHERE\n" +
                        "    (position >= %d AND position <= %d) AND\n" +
                        "    position IN\n" +
                        "    (\n" +
                        "      SELECT DISTINCT position FROM `%s`\n" +
                        "      WHERE position >= %d AND position <= %d\n" +
                        "    )\n" +
                        ")\n" +
                        "SELECT\n" +
                        "  new_pet.position,\n" +
                        "  ARRAY_AGG(STRUCT(\n" +
                        "    new_pet.sample,\n" +
                        "    state,\n" +
                        "    ref,\n" +
                        "    alt,\n" +
                        "    AS_RAW_MQ,\n" +
                        "    AS_RAW_MQRankSum,\n" +
                        "    AS_QUALapprox,\n" +
                        "    AS_RAW_ReadPosRankSum,\n" +
                        "    AS_SB_TABLE,\n" +
                        "    AS_VarDP,\n" +
                        "    call_GT,\n" +
                        "    call_AD,\n" +
                        "    call_DP,\n" +
                        "    call_GQ,\n" +
                        "    call_PGT,\n" +
                        "    call_PID,\n" +
                        "    call_PL\n" +
                        "  )) AS values\n" +
                        "FROM\n" +
                        "  new_pet\n" +
                        "LEFT OUTER JOIN\n" +
                        "  (SELECT * from `%s` AS vet_inner\n" +
                        "   WHERE vet_inner.position >= %d AND vet_inner.position <= %d)\n" +
                        "  AS vet\n" +
                        "ON (new_pet.position = vet.position AND new_pet.sample = vet.sample)\n" +
                        "GROUP BY new_pet.position\n" +
                        limitString,

                getFQPositionTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd()
        );
    }

    private String getOptimizedUngroupedVariantQueryString(final SimpleInterval interval) {
        String limitString = "";
        if (queryRecordLimit > 0) {
            limitString = "LIMIT " + queryRecordLimit;
        }

        return String.format(
                "WITH new_pet AS\n" +
                        "(\n" +
                        "  SELECT * FROM `%s`\n" +
                        "  WHERE\n" +
                        "    (position >= %d AND position <= %d) AND\n" +
                        "    position IN\n" +
                        "    (\n" +
                        "      SELECT DISTINCT position FROM `%s`\n" +
                        "      WHERE position >= %d AND position <= %d\n" +
                        "      AND sample IN (SELECT sample FROM `%s`)\n" +
                        "    )\n" +
                        "    AND sample IN (SELECT sample FROM `%s`)\n" +
                        ")\n" +
                        "SELECT\n" +
                        "  new_pet.position,\n" +
                        "  new_pet.sample,\n" +
                        "  state,\n" +
                        "  ref,\n" +
                        "  alt,\n" +
                        "  AS_RAW_MQ,\n" +
                        "  AS_RAW_MQRankSum,\n" +
                        "  AS_QUALapprox,\n" +
                        "  AS_RAW_ReadPosRankSum,\n" +
                        "  AS_SB_TABLE,\n" +
                        "  AS_VarDP,\n" +
                        "  call_GT,\n" +
                        "  call_AD,\n" +
                        "  call_DP,\n" +
                        "  call_GQ,\n" +
                        "  call_PGT,\n" +
                        "  call_PID,\n" +
                        "  call_PL\n" +
                        "FROM\n" +
                        "  new_pet\n" +
                        "LEFT OUTER JOIN\n" +
                        "  (SELECT * from `%s` AS vet_inner\n" +
                        "   WHERE vet_inner.position >= %d AND vet_inner.position <= %d)\n" +
                        "  AS vet\n" +
                        "ON (new_pet.position = vet.position AND new_pet.sample = vet.sample)\n" +
                        limitString,

                getFQPositionTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQTableName(contigToSampleTableMap.entrySet().iterator().next().getValue()),
                getFQTableName(contigToSampleTableMap.entrySet().iterator().next().getValue()),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd()
        );
    }

    private String getCohortExtractQueryString(final SimpleInterval interval) {
        String limitString = "";
        if (queryRecordLimit > 0) {
            limitString = "LIMIT " + queryRecordLimit;
        }

        return String.format(
                "WITH new_pet AS\n" +
                        "(\n" +
                        "  SELECT * FROM `%s`\n" +
                        "  WHERE\n" +
                        "    (position >= %d AND position <= %d) AND\n" +
                        "    position IN\n" +
                        "    (\n" +
                        "      SELECT DISTINCT position FROM `%s`\n" +
                        "      WHERE position >= %d AND position <= %d\n" +
                        "      AND sample IN (SELECT sample FROM `%s`)\n" +
                        "    )\n" +
                        "    AND sample IN (SELECT sample FROM `%s`)\n" +
                        ")\n" +
                        "SELECT\n" +
                        "  new_pet.position,\n" +
                        "  new_pet.sample,\n" +
                        "  new_pet.state,\n" +
                        "  vet.ref,\n" +
                        "  REPLACE(vet.alt,\",<NON_REF>\",\"\") alt,\n" +
                        "  yng.vqslod  call_AS_VQSLOD,\n" +
                        "  yng.yng_status call_YNG_STATUS,\n" +
                        "  call_GT,\n" +
                        "  call_GQ,\n" +
                        "  cast(SPLIT(call_pl,\",\")[OFFSET(0)] as int64) as call_RGQ\n" +
                        "FROM\n" +
                        "  new_pet\n" +
                        "LEFT OUTER JOIN\n" +
                        "  (SELECT * from `%s` AS vet_inner\n" +
                        "   WHERE vet_inner.position >= %d AND vet_inner.position <= %d)\n" +
                        "  AS vet\n" +
                        "ON (new_pet.position = vet.position AND new_pet.sample = vet.sample)\n" +
                        "LEFT OUTER JOIN\n" +
                        "  `%s` AS yng\n" +
                        "ON (new_pet.position = yng.position AND vet.ref = yng.ref AND REPLACE(vet.alt,\",<NON_REF>\",\"\") = yng.alt)" +
                        limitString,

                getFQPositionTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQTableName(contigToSampleTableMap.entrySet().iterator().next().getValue()),
                getFQTableName(contigToSampleTableMap.entrySet().iterator().next().getValue()),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd(),
                filteringFQTableName
        );
    }

    private String getVQSRFeatureExtractQueryString(final SimpleInterval interval, final boolean trainingSitesOnly) {
            String trainingSitesStanza =
                    !trainingSitesOnly?"":
                            "AND position IN (SELECT position FROM `broad-dsp-spec-ops.joint_genotyping_ref.vqsr_training_sites_*` WHERE chrom='chr20')\n";
            String query =
                    "WITH ref_ad_info AS (\n" +
                            "SELECT \n" +
                            "  position,\n" +
                            "  IFNULL(sum(cast(SPLIT(call_AD,\",\")[OFFSET(0)] as int64)),0) as ref_ad\n" +
                            "FROM `@vet` \n" +
                            "WHERE sample IN (SELECT sample FROM `@sample`)\n" +
                            "AND (position >= @start AND position <= @end) \n" +
                            trainingSitesStanza +
                            "AND (\n" +
                            "  SELECT SUM(CAST(part AS int64)) FROM UNNEST(SPLIT(call_AD, ',')) part WITH OFFSET index WHERE index >= 1 ) \n" +
                            "  > 1\n" +
                            "GROUP BY position),\n" +
                            "ref_sb_info AS (\n" +
                            "SELECT \n" +
                            "  position,\n" +
                            "  IFNULL(sum(cast(SPLIT(SPLIT(as_sb_table,\"|\")[OFFSET(0)],\",\")[OFFSET(0)] as int64)),0) as sb_ref_plus, \n" +
                            "  IFNULL(sum(cast(SPLIT(SPLIT(as_sb_table,\"|\")[OFFSET(0)],\",\")[OFFSET(1)] as int64)),0) as sb_ref_minus  \n" +
                            "FROM `@vet` \n" +
                            "WHERE sample IN (SELECT sample FROM `@sample`)\n" +
                            "AND (position >= @start AND position <= @end) \n" +
                            trainingSitesStanza +
                            "GROUP BY position)\n" +
                            "SELECT \n" +
                            "       ai.position, \n" +
                            "       ai.ref, \n" +
                            "       ai.allele, \n" +
                            "       RAW_QUAL,\n" +
                            "       radi.ref_ad,\n" +
                            "       AS_MQRankSum,\n" +
                            "       AS_MQRankSum_ft,\n" +
                            "       AS_ReadPosRankSum,\n" +
                            "       AS_ReadPosRankSum_ft,\n" +
                            "       RAW_MQ,\n" +
                            "       RAW_AD,\n" +
                            "       RAW_AD_GT_1,\n" +
                            "       rsbi.SB_REF_PLUS,\n" +
                            "       rsbi.SB_REF_MINUS,\n" +
                            "       SB_ALT_PLUS, \n" +
                            "       SB_ALT_MINUS\n" +
                            "FROM (\n" +
                            "SELECT aa.position, \n" +
                            "       ref, \n" +
                            "       allele, \n" +
                            "       IFNULL(SUM(qual),0) as RAW_QUAL,\n" +
                            "       `bqutil`.fn.median(ARRAY_AGG( raw_mqranksum_x_10 IGNORE NULLS)) / 10.0 as AS_MQRankSum,\n" +
                            "       `broad-dsp-spec-ops`.joint_genotyping_ref.freq_table(ARRAY_AGG(raw_mqranksum_x_10 IGNORE NULLS)) AS_MQRankSum_ft,\n" +
                            "       `bqutil`.fn.median(ARRAY_AGG(raw_readposranksum_x_10 IGNORE NULLS)) / 10.0 as AS_ReadPosRankSum,\n" +
                            "       `broad-dsp-spec-ops`.joint_genotyping_ref.freq_table(ARRAY_AGG(raw_readposranksum_x_10 IGNORE NULLS)) as AS_ReadPosRankSum_ft,\n" +
                            "       IFNULL(SUM(RAW_MQ),0) as RAW_MQ,\n" +
                            "       IFNULL(SUM(AD),0) as RAW_AD, \n" +
                            "       IFNULL(SUM(CASE WHEN AD > 1 THEN AD ELSE 0 END),0) as RAW_AD_GT_1, # to match GATK implementation\n" +
                            "       IFNULL(SUM(SB_ALT_PLUS),0)  as SB_ALT_PLUS, \n" +
                            "       IFNULL(SUM(SB_ALT_MINUS),0) as SB_ALT_MINUS\n" +
                            "FROM `@altAllele` as aa\n" +
                            "WHERE (position >= @start AND position <= @end) \n" +
                            trainingSitesStanza +
                            "AND sample in (SELECT sample from `@sample` )\n" +
                            "AND allele != '*'\n" +
                            "GROUP BY 1,2,3\n" +
                            ") ai\n" +
                            "LEFT JOIN ref_ad_info radi ON (ai.position = radi.position)\n" +
                            "LEFT JOIN ref_sb_info rsbi ON (ai.position = rsbi.position)\n";

            return query
                .replaceAll("@vet", getFQVariantTable(interval))
                .replaceAll("@pet", getFQPositionTable(interval))
                .replaceAll("@sample", getFQTableName(contigToSampleTableMap.entrySet().iterator().next().getValue()))
                .replaceAll("@start", String.format("%d",interval.getStart()))
                .replaceAll( "@end", String.format("%d",interval.getEnd()))
                .replaceAll( "@altAllele", getFQTableName(contigToAltAlleleTableMap.entrySet().iterator().next().getValue()));
    }

    private String getSampleListQueryString(final String sampleTableName) {
        return "SELECT sample FROM `" + getFQTableName(sampleTableName)+ "`";
    }
}
